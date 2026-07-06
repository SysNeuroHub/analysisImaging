%% this is a script to do the following:
% load data from Market server
% detect onset of each stimulation in thorsync time using photodiode signal
% detect wf frame time in thorsync time
% obtain stimulus triggered time courses
% compute Direction selectivity index of each cell



setPath_analysisImaging;
saveDir = '/home/daisuke/Documents/git/analysisImaging/wfAnalysis/hctsa_tmp';

%% experiment

expt.subject = 'Confucious';
expt.expDate = '2026-04-25_1';
expt.expNum = 6;
bklightCtrl = 0;

%% SVD
nSV = 1000;

%% analysis
marginT = .2; %[s]
resizeS = 0.5; %spatial rescaling factor

params.movieSuffix = 'amber';
params.useCorrected = 0;


%% data location
thisDate = expt.expDate(1:10);
thisSeries = str2num(expt.expDate(12:end));
expPath = fileparts(dat.expFilePath(expt.subject, thisDate, thisSeries, expt.expNum, 'widefield','master'));
saveVpath = expPath;
mpepDir = dat.reposPath('main', 'master');

load(dat.expFilePath(expt.subject, thisDate, thisSeries, expt.expNum, 'Timeline', 'master'));


%% load wf data

tname = params.movieSuffix;
if params.useCorrected
    tname = [tname '-corrected'];
end

disp('Loading widefield data');
disp(expt)
[U, V, t, mimg, mask] = quickLoadUVt(expPath, nSV, saveVpath, params);
if isempty(t) %HACK when timeline strobes did not align with tiff frames
    [strobeOnTimes, strobeOffTimes, strobeDurs] = getStrobeTimes(Timeline, 'alloptrig');
    t = strobeOnTimes;
    t = t(1:size(V,2));
end

Fs = 1/median(diff(t));

U = imresize(U,resizeS);
mask = imresize(mask, resizeS);
mimg = imresize(mimg, resizeS);
Ux = size(U,2);
Uy = size(U,1);

%% temporal filtering V
fV = filtV(V, Fs, 0.01, 5);

%% load mpep data
p = ProtocolLoad_wf(expt.subject,expt.expDate,expt.expNum); %3/6/20
figname = [dat.constructExpRef(expt.subject, thisDate, thisSeries, expt.expNum) '_'...
    p.xfile(1:end-2) '_' params.movieSuffix '_' num2str(params.useCorrected)];


if contains(p.xfile, 'stimORsequence')
    doSequence = 1;
else
    doSequence = 0;
end

stimSequence = getStimSequence(expt.subject, expt.expDate, expt.expNum);
%stimSequence.seq
%stimSequence.labels
%stimSequence.paramValues
if ~isfield(p,'pfiledurs')
    p.pfiledurs = p.pars(1,:)/10;
end

saveName = [dat.constructExpRef(expt.subject, thisDate, thisSeries, expt.expNum) '_'...
        p.xfile(1:end-2) '_' params.movieSuffix '_' num2str(params.useCorrected) '_singleTraces'];

respWin = [-marginT min(p.pfiledurs)+marginT]; %31/3/20


%% detect onset of each stimulation in Timeline time using photodiode signal
expt = grabStimTimesWF(expt, 0, [], [], bklightCtrl,0);
laserTh = .5; %[V]
expt = grabLaserTimesWF(expt,[],laserTh);



%check number of stimulation
if length(expt.stimTimes.onset) ~= p.nstim * p.nrepeats
    error(['Detected #stim onset:' num2str(expt.stimTimes.onset') '/ total #stim' num2str(p.nstim * p.nrepeats)]);
end

%% stimulus triggered movie
% pixelTuningCurveViewerSVD(U, V, t, expt.laserTimes.onset, stimSequence.seq, respWin, 1);
% title(tname);


dur = 10;%[s]
nEpochs = 20;
mrec = [];
for iepoch = 1:nEpochs
    trange = [dur*iepoch dur*(iepoch+1)];
    [~,tidx(1)] = min(abs(t-trange(1)));
    [~,tidx(2)] = min(abs(t-trange(2)));
    rec = svdFrameReconstruct(U,V(:,tidx(1):tidx(2)));
    mrec(:,iepoch) = squeeze(nanmean(nanmean(rec,1),2));
end


save(fullfile(saveDir,saveName), 'dur','mrec');
