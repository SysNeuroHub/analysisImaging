%% this is a script to do the following:
% load data from Market server
% detect onset of each stimulation in thorsync time using photodiode signal
% detect wf frame time in thorsync time
% obtain stimulus triggered time courses
% compute Direction selectivity index of each cell



setPath_analysisImaging;

%% experiment

expt.subject = 'test_stimTTLOsc';
expt.expDate = '2026-02-13_2';
expt.expNum = 1;
bklightCtrl = 0;

%% SVD
nSV = 1000;

%% analysis
marginT = .2; %[s]
resizeS = 0.5; %spatial rescaling factor

params.movieSuffix = 'amber';
params.useCorrected = 0;

resultSaveDir = fullfile('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\AnalysisResult\wf',...
    [expt.subject '_' expt.expDate '_' num2str(expt.expNum)]);
mkdir(resultSaveDir);


%% data location
thisDate = expt.expDate(1:10);
thisSeries = str2num(expt.expDate(12:end));
expPath = fileparts(dat.expFilePath(expt.subject, thisDate, thisSeries, expt.expNum, 'widefield','master'));
saveVpath = expPath;
mpepDir = dat.reposPath('main', 'master');

load(dat.expFilePath(expt.subject, thisDate, thisSeries, expt.expNum, 'Timeline', 'master'));


%% load wf data
% for icolor = 2
%
%     switch icolor
%         case 1
%             params.movieSuffix = 'red';
%             params.useCorrected = 0;
%         case 2
%             params.movieSuffix = 'amber';
%             params.useCorrected = 0;
%         case 3
%             params.movieSuffix = 'amber';
%             params.useCorrected = 1;
%     end
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

%% dF/F NG
[newU, newV] = dffFromSVD(U, V, mimg);

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
pixelTuningCurveViewerSVD(U, V, t, expt.stimTimes.onset, stimSequence.seq, respWin, 0);
title(tname);



%% frequency response
ComplexMap = []; AbsMap = []; AngleMap = [];
for icond = 1:p.nstim
    these = find(stimSequence.seq == icond);
    stimFreq = 2*p.pars(5,icond)/10; %TEMP
    stimOn = p.pars(3,icond)*1e-3;
    stimOff = p.pars(4,icond)*1e-3;
    
    [avgPeriEventV, winSamps, periEventV, sortLabels] = ...
        eventLockedAvg(V, t, expt.stimTimes.onset(these), ...
        stimSequence.seq(these), [stimOn stimOff]);
    
    [ComplexMap(:,:,icond), AbsMap(:,:,icond), AngleMap(:,:,icond)] = ...
        FourierMapSVD(U, squeeze(avgPeriEventV), winSamps, stimFreq);
end

images(AbsMaps);

%% single trial traces
%     yy=189; xx = 177;
%     icond = 17;
%     these = find(sortLabels == icond);
%     trace = [];
%     for ii = 1:numel(these)
%         trace(:,ii)=svdFrameReconstruct(U(yy,xx,:),squeeze(periEventV(these(ii),:,:)));
%     end
%     imagesc(winSamps, 1:numel(these), trace');


