%% this is a script to do the following:
% load data from Market server
% detect onset of each stimulation in thorsync time using photodiode signal
% detect wf frame time in thorsync time
% obtain stimulus triggered time courses
% compute Direction selectivity index of each cell



setPath_analysisImaging;

%% experiment
expt.subject = 'yamatotakeru';
expt.expDate = '2025-04-21_1';
expt.expNum = 12;
bklightCtrl = 0;

%% SVD
nSV = 1000;%1000;
params.movieSuffix = 'amber';% 'purple'
params.useCorrected = 1;

%% analysis
marginT = .5; %[s]
resizeS = 0.25; %spatial rescaling factor



%for estimation of preferred stim
n_boot = 10;%1 to see retinotopy only, 10 to compute VFS
use_method = 'max'; % max or com
screen_resize_scale = 3; %3 if max method


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
fV = filtV(V, Fs, 0.01, 2);

%% load mpep data
p = ProtocolLoad_wf(expt.subject,expt.expDate,expt.expNum); %3/6/20
figname = [dat.constructExpRef(expt.subject, thisDate, thisSeries, expt.expNum) '_'...
    p.xfile(1:end-2) '_' params.movieSuffix];


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
expt = grabStimTimesWF(expt, 0, [], [], bklightCtrl,1);
% expt.stimTimes.onset
% expt.stimTimes.offset
% expt.stimTimes.frameTimes
%laserTh = 0.1;
%expt = grabLaserTimesWF(expt,[],[],[],laserTh);



%check number of stimulation
if length(expt.stimTimes.onset) ~= p.nstim * p.nrepeats
    error(['Detected #stim onset:' num2str(expt.stimTimes.onset') '/ total #stim' num2str(p.nstim * p.nrepeats)]);
end

% %% hack for izanami 2025-01-30_2 exp4
% lastValidTime = expt.stimTimes.offset(18*3);
% [~, lastValidTidx] = min(abs(t - lastValidTime));
% 
% t = t(1:lastValidTidx);
% V = V(:,1:lastValidTidx);
% fV = fV(:,1:lastValidTidx);
% newV = newV(:,1:lastValidTidx);

%% stimulus triggered movie
pixelTuningCurveViewerSVD(U, V, t, expt.stimTimes.onset, stimSequence.seq, respWin);
[avgPeriEventV, winSamps, periEventV] = ...
    eventLockedAvg(V, t, expt.stimTimes.onset, stimSequence.seq, respWin);
%avgPeriEventV: icond x nSV x time
%periEventV: event x nSV x time



    %% time-avg response & stimulus preference map
    preIdx = find(winSamps<0);
     [~,postIdx] = min(abs(winSamps - 0.25));%to examine laser artifact
%      postIdx = intersect(find(winSamps>1), find(winSamps < 2));%min(p.pfiledurs)));

%subtract by prestimulus in each condition
tavgRespV = mean(avgPeriEventV(:,:,postIdx),3) - mean(avgPeriEventV(:,:,preIdx),3);
%condition x nSV

tavgResp = svdFrameReconstruct(U, tavgRespV');
%tavgResp = tavgResp - tavgResp(:,:,p.blankstims);%still blood vessel remains...

nRows            = 1;%ceil(sqrt(p.nstim));
nCols = p.nstim-1;%ceil(p.nstim/nRows);
figure('position',[0 0 1900 1200]);
panel = [];
for istim = 1:p.nstim-1
    panel(istim) = subplot(nRows,nCols,istim);
    imagesc(100*squeeze(tavgResp(:,:,istim)./mimg),'alphadata',mask);
    axis equal tight;
    title(stimSequence.labels{istim});
   caxis([-3 3])
    [h,g]=mcolorbar;
end
g.YLabel.String='dF/F [%]';
saveas(gcf,fullfile(resultSaveDir,figname),'fig');
close;


%% show movie
traces = prepareTimelineTraces(Timeline);
movieWithTracesSVD(U, fV, t, traces);

