%% this is a script to do the following:
% load data from Market server
% detect onset of each stimulation in thorsync time using photodiode signal
% detect wf frame time in thorsync time
% obtain stimulus triggered time courses
% compute Direction selectivity index of each cell



setPath_analysisImaging;

%% experiment
expt.subject = 'robita';
expt.expDate = '2025-03-29_1';
expt.expNum = 1;%2;
bklightCtrl = 0;

%% SVD
nSV = 1000;%1000;
params.movieSuffix = 'amber';% 'purple'
params.useCorrected = 0;

%% analysis
marginT = .5; %[s]
resizeS = 0.25; %spatial rescaling factor



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

%% stimulus triggered movie
pixelTuningCurveViewerSVD(U, V, t, expt.stimTimes.onset, stimSequence.seq, respWin);
[avgPeriEventV, winSamps, periEventV] = ...
    eventLockedAvg(V, t, expt.stimTimes.onset, stimSequence.seq, respWin);
%avgPeriEventV: icond x nSV x time
%periEventV: event x nSV x time


istim = 5;
t1 = 0;%2.5;
t2 = 2.5;%3;
ypixRange = [181 220];  
xpixRange = [51 90];
crange = [-.054 .054];

%% time-avg response & stimulus preference map
preIdx = find(winSamps<0);
postIdx = intersect(find(winSamps>t1), find(winSamps < t2));

%dF/F
avgResp = 100*svdFrameReconstruct(U, squeeze(avgPeriEventV(istim,:,:)))./repmat(mimg,[1 1 numel(winSamps)]);
avgResp = avgResp - mean(avgResp(:,:,preIdx),3);
tavgResp = mean(avgResp(:,:,postIdx),3);

figure('position',[0 0 1000 300]);
panel = [];
panel(istim) = subplot(121);
imagesc(tavgResp,'alphadata',mask);
axis equal tight off;
hold on;
rectangle('position', [xpixRange(1) ypixRange(1) diff(xpixRange) diff(ypixRange)], 'curvature',1);
%title(stimSequence.labels{istim});
title([num2str(t1) '-' num2str(t2) 's after optoStim onset']);
caxis(crange)
colormap(flipud(parula));
[h,g]=mcolorbar(gca,.5);
g.YLabel.String='dF/F [%]';

MmPerPixel_t = 0.0104 / resizeS;

% bregma = [110 135];
% lambda = [227 135];
% addAllenCtxOutlines(bregma, lambda, 'w', MmPerPixel_t);
%drawTopDownCtx

subplot(122);
plot(winSamps, squeeze(mean(mean(avgResp(ypixRange(1):ypixRange(2),xpixRange(1):xpixRange(2),:)))));
axis ij; ylim(crange);
set(gca,'tickdir','out')
vbox([0 1 2],[.5 1.5 2.5],gca,[.8 .8 1]);
vbox(t1,t2,gca,[1 .8 .8]);hline(0);
xlabel('Time [s]'); ylabel('dF/F[%]');
xlim([-0.5 8.5]);

savePaperFigure(gcf,[figname '_icond' num2str(istim)]);


