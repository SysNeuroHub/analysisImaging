% created from retinotopyKalatsky_singleTr_bootstrap.m 10/1/2016

addpath(genpath('C:\npy-matlab'));
addpath(genpath('C:\Users\dshi0006\npy-matlab'));
addpath(genpath('C:\Users\Analysis\npy-matlab'));
addpath(genpath('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\master\Analysis\wfAnalysis\widefield'));
addpath(genpath('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\master\visbox\rigbox'));
addpath('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\master\Analysis');
addpath('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\master\Analysis\wfAnalysis');


%% experiment
%ppbox notation
expt.subject = 'L4GCaMP6s_260';
expt.expDate = '2020-10-14_1';
expt.expNum = 2;
bklightCtrl = 0;

%% SVD
nSV = 100;
params.movieSuffix = 'blue';% 'purple''corr_dFF'; %'blue'  %U for corr_dFF can have NANs..
params.useCorrected = 0;

%% analysis
resizeS = 0.25; %spatial rescaling factor
stimList = [1:4];
filterFourierMaps = 0;
filtRF = 2;%spatial filter size in pixels

resultSaveDir = fullfile('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\AnalysisResult\wf',...
    [expt.subject '_' expt.expDate '_' num2str(expt.expNum)]);
mkdir(resultSaveDir);


%% data location
thisDate = expt.expDate(1:10);
thisSeries = str2num(expt.expDate(12:end));
expPath = fileparts(dat.expFilePath(expt.subject, thisDate, thisSeries, expt.expNum, 'widefield','master'));
saveVpath = expPath;
mpepDir = dat.reposPath('main', 'master');


%% load wf data
disp('Loading widefield data');
disp(expt)
[U, V, t, mimg, mask] = quickLoadUVt(expPath, nSV, saveVpath, params);
Fs = 1/median(diff(t));


U = imresize(U,resizeS);
mask = imresize(mask, resizeS);
Ux = size(U,2);
Uy = size(U,1);

%movieWithTracesSVD(U,V,t)
%svdViewer(U, Sv, V, 35);

% %% single trial traces
% % TODO: make a GUI to pick a pixel of interest
% roiy = 73:74;%320:325;%round(Ux/2)-2:round(Ux/2)+2;
% roix = 118:119;%100:105;%round(Uy/2)-2:round(Uy/2)+2;
% 
% % TODO:
% % disp('select pixel');
% % [roix, roiy] = pickPixel(mimg);
% % keyboard;
% 
% ROItrace = squeeze(mean(mean(svdFrameReconstruct(U(roiy,roix,:),V))));
% irepeat = 10;
% 
% figure('position',[0 0 1000 1000]);
% showSingleTrialTraces(ROItrace, irepeat, t, expt, p, stimSequence, respWin, prctile(ROItrace,[1 99.8]));
% %%

%% load screenInfo
expt = grabScreenInfo(expt);

%% load p-file
p = ProtocolLoad_wf(expt.subject,expt.expDate,expt.expNum); %3/6/20
stimSequence = getStimSequence(expt.subject, expt.expDate, expt.expNum);
%stimSequence.seq
%stimSequence.labels
%stimSequence.paramValues


%% load parameters on stimulus
pval_ori = p.pars(find(strcmp(p.parnames, 'ori')),:);
pval_maxL = p.pars(find(strcmp(p.parnames, 'maxL')),:); %empty
stimidxs_h = intersect(find(pval_ori == 1), find(pval_maxL ~= 0)); %empty
stimidxs_v = intersect(find(pval_ori == 2), find(pval_maxL ~= 0)); %empty

pidx_start = find(strcmp(p.parnames, 'start'));
pidx_end = find(strcmp(p.parnames, 'end'));

start_point_h = min([p.pars(pidx_start,stimidxs_h) p.pars(pidx_end,stimidxs_h)]);%empty
end_point_h = max([p.pars(pidx_start,stimidxs_h) p.pars(pidx_end,stimidxs_h)]);%empty
start_point_v = min([p.pars(pidx_start,stimidxs_v) p.pars(pidx_end,stimidxs_v)]);%empty
end_point_v = max([p.pars(pidx_start,stimidxs_v) p.pars(pidx_end,stimidxs_v)]);%empty

pidx_tf = find(strcmp(p.parnames, 'tf'));
tf = unique(p.pars(pidx_tf,stimList)) / 100; %[Hz]

nFramesPerSweep = round(expt.screenInfo.FrameRate/tf);
tf_corrected = expt.screenInfo.FrameRate / nFramesPerSweep;

%% get stimulus onset/offset from Timeline
expt = grabStimTimesWF(expt, 0,[],[],bklightCtrl);

pixelTuningCurveViewerSVD(U, V, t, expt.stimTimes.onset, stimSequence.seq, [0 200], 0);

%% frequency analysis
%get phase and amplitude at tf_corrected

ANGLEMAPS = zeros(Uy,Ux,4);
ABSMAPS = zeros(Uy,Ux,4);
for iStim = 1:4
    
    [periEventV, winSamps] = eventLockedAvg(V, t, expt.stimTimes.onset(iStim), ...
        stimSequence.seq(iStim), [0 expt.stimTimes.offset(iStim)-expt.stimTimes.onset(iStim)]);
    periEventV = squeeze(periEventV);%omit the 1st dimension
    
    
    Tensor =  svdFrameReconstruct(U, periEventV);
    %should subtract by mimg or else?
    %do detrend?
    
    [ComplexMaps, AbsMaps, AngleMaps] = FourierMaps(Tensor, Fs, tf_corrected);
    
    if filterFourierMaps
        realMaps_filt = wiener2(real(ComplexMaps), [filtRF filtRF]);
        imagMaps_filt = wiener2(imag(ComplexMaps), [filtRF filtRF]);
        
        ComplexMaps_filt = realMaps_filt + 1i * imagMaps_filt;
        AngleMaps = angle(ComplexMaps_filt);
    end
    %         %correction of angle produced by stackset.fouriermaps. range = [0 2pi]
    %         AngleMaps = -AngleMaps_filt;
    %         AngleMaps(AngleMaps>2*pi) = AngleMaps(AngleMaps>2*pi) - 2*pi;
    %         AngleMaps(AngleMaps<startPhase(iStim)) = AngleMaps(AngleMaps<startPhase(iStim)) + 2*pi;
    
    ANGLEMAPS(:, :,  stimSequence.seq(iStim)) = AngleMaps;
    ABSMAPS(:,:,stimSequence.seq(iStim)) = AbsMaps;
end


%% combine opposite directions ... needs revision
xMap = ANGLEMAPS(:,:,1) - ANGLEMAPS(:,:,2); %needs revision 
yMap = ANGLEMAPS(:,:,3) - ANGLEMAPS(:,:,4); %needs revision


%% compute visual field sign map from xMap and yMap
signMap = VFS(xMap, yMap);



%% visualize
figname = [dat.constructExpRef(expt.subject, thisDate, thisSeries, expt.expNum) '_'...
    p.xfile(1:end-2) '_' params.movieSuffix];

for iStim = 1:4
    subplot(2,2,iStim);
    imagesc(ANGLEMAPS(:,:,iStim));
    caxis([-pi pi]);
    colormap(hsv);
    title(stimSequence.labels{iStim});
end
mcolorbar(gcf, .5);
screen2png(fullfile(resultSaveDir,[figname '_' params.movieSuffix 'angleMaps.png']));
close;


for iStim = 1:4
    subplot(2,2,iStim);
    imagesc(ABSMAPS(:,:,iStim));
    title(stimSequence.labels{iStim});
end
mcolorbar(gcf, .5);
screen2png(fullfile(resultSaveDir,[figname '_' params.movieSuffix 'absMaps.png']));
close;

%todo: single trial response at specified pixel as in mainScript_wf.m

subplot(131);imagesc(xMap);axis equal tight; 
caxis(prctile(xMap(:),[10 90]));mcolorbar(gca,.5);
grid on

subplot(132);imagesc(yMap);axis equal tight; 
caxis(prctile(yMap(:),[10 90]));mcolorbar(gca,.5);
grid on

subplot(133);imagesc(signMap);axis equal tight; 
caxis([-1 1]);mcolorbar(gca,.5);
grid on

screen2png(fullfile(resultSaveDir,[figname '_' params.movieSuffix 'signMap.png']));
close;