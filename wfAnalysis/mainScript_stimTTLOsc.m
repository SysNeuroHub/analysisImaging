%% this is a script to do the following:
% load data from Market server
% detect onset of each stimulation in thorsync time using photodiode signal
% detect wf frame time in thorsync time
% obtain stimulus triggered time courses
% compute Direction selectivity index of each cell



setPath_analysisImaging;

%% experiment

expt.subject = 'Confucious';
expt.expDate = '2026-03-24_1';
expt.expNum = 2;
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
pixelTuningCurveViewerSVD(U, V, t, expt.laserTimes.onset, stimSequence.seq, respWin, 1);
title(tname);



%% frequency response
ComplexMap = []; AbsMap = []; AngleMap = [];
for icond = 1:p.nstim
    these = find(stimSequence.seq == icond);
    stimFreq = 2*p.pars(5,icond)/1000;

    for iphase = 1:3
        switch iphase
            case 1
                On = 0;
                Off = p.pars(3,icond)*1e-3;
            case 2
                On = p.pars(3,icond)*1e-3;
                Off = p.pars(4,icond)*1e-3;
            case 3
                On = p.pars(4,icond)*1e-3;
                Off = p.pars(1,icond)*1e-1;
        end
        [avgPeriEventV, winSamps, periEventV, sortLabels] = ...
            eventLockedAvg(V, t, expt.stimTimes.onset(these), ...
            stimSequence.seq(these), [On Off]);

        [ComplexMap(:,:,icond,iphase), AbsMap(:,:,icond,iphase), AngleMap(:,:,icond,iphase)] = ...
            FourierMapSVD(U, squeeze(avgPeriEventV), winSamps, stimFreq);
    end

    subplot(p.nstim, 2, 2*icond-1);imagesc(AbsMap(:,:,icond,2)./AbsMap(:,:,icond,1)); clim([0 2]); axis equal tight 
    if icond==1; title('amplitude during stimulation/prestim'); end
    ylabel(['icond' num2str(icond) ', ' num2str(stimFreq) ' Hz']);
    subplot(p.nstim, 2, 2*icond);imagesc(AbsMap(:,:,icond,3)./AbsMap(:,:,icond,1)); clim([0 2]); axis equal tight; mcolorbar
    if icond==1; title('amplitude during poststim/prestim'); end
end

%% single trial traces
icond = 3;
these = find(stimSequence.seq == icond);
stimFreq = p.pars(5,icond)/1000;
dur = p.pars(1,icond)*1e-1;
stimOn = p.pars(3,icond)*1e-3;
stimOff = p.pars(4,icond)*1e-3;

[avgPeriEventV, winSamps, periEventV, sortLabels] = ...
    eventLockedAvg(V, t, expt.stimTimes.onset(these), ...
    stimSequence.seq(these), [0 dur]);

yy=171; xx = 131;
trace = [];
for ii = 1:numel(these)
    trace(:,ii)=svdFrameReconstruct(U(yy,xx,:),squeeze(periEventV(ii,:,:)));
end
imagesc(winSamps, 1:numel(these), trace');
xlabel('time [s]');
ylabel('trial ID');
title(['yy=' num2str(yy) ', xx=' num2str(xx) ', stimulus frequency ' num2str(stimFreq)  ' [Hz]']);
vline([stimOn stimOff]);
mcolorbar;


