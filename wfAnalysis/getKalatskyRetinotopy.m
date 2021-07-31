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
expt.subject = 'L4GCaMP6s_252';
expt.expDate = '2021-01-07_1';
expt.expNum = 2;
bklightCtrl = 0;

%% SVD
nSV = 1000;
params.movieSuffix = 'blue';% 'purple''corr_dFF'; %'blue'  %U for corr_dFF can have NANs..
params.useCorrected = 1;

%% analysis
resizeS = 0.5; %spatial rescaling factor
stimList = [1:4];
filterFourierMaps = 1;
filtRF = 3;
delayPhase = 0;%-pi/2;%-3*pi/4;%pi; %added delay in radian
%0 for intrinsic imaging
%pi for calcium imaging
n_boot = 10;

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

% %% align CCF
% expt = grabStereoInfo(expt);
% expt.stereoInfo.bregma = resizeS*(expt.stereoInfo.bregma);
% expt.stereoInfo.lambda = resizeS*(expt.stereoInfo.lambda);
% 
% %TODO: retrieve this from ops.mat?
% objectiveType = 0.5;
% hwbinning = 2;
% expt.mmPerPixel = 1e-3*5/objectiveType*hwbinning/resizeS;
% 
% imagesc(mimg);
% addAllenCtxOutlines(expt.stereoInfo.bregma, expt.stereoInfo.lambda,'r',expt.mmPerPixel);


%% temporal filtering 14/10/20 ... cause some weird artefact on restored trace
% highpassCutoff = 0.01;
% heartbeatBandStop = [3 5];
% V = detrendAndFilt(V, Fs, highpassCutoff, heartbeatBandStop); 


U = imresize(U,resizeS);
mimg = imresize(mimg, resizeS);
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

if strcmp(expt.subject, 'TIGRE2GCaMP6s_318') && strcmp(expt.expDate, '2020-11-05_3') ...
        && isequal(expt.expDate,4)
    p.nrepeats=3; %TEMP
end

%% load parameters on stimulus
pval_ori = p.pars(find(strcmp(p.parnames, 'ori')),:);

pidx_tf = find(strcmp(p.parnames, 'tf'));
tf = unique(p.pars(pidx_tf,stimList)) / 100; %[Hz]

nFramesPerSweep = round(expt.screenInfo.FrameRate/tf);
tf_corrected = expt.screenInfo.FrameRate / nFramesPerSweep;

dur = unique(p.pars(find(strcmp(p.parnames, 'dur'))))/10;%[s]
blankDur = unique(p.pars(find(strcmp(p.parnames, 'bdur'))))/1000;%[s]
stimDur = 1/tf_corrected-blankDur; %[s] %duration when stimulus is presented in one cycle

%% get stimulus onset/offset from Timeline
ISI = 0.8*p.minWait; %15/10/20 hack
expt = grabStimTimesWF(expt, 0, [], ISI, bklightCtrl, 0);

% pixelTuningCurveViewerSVD(U, V, t, expt.stimTimes.onset, stimSequence.seq, [-5 dur]);, 1);

%% frequency analysis
%get phase and amplitude at tf_corrected
[avgPeriEventV, winSamps, periEventV] = eventLockedAvg(V, t, expt.stimTimes.onset, ...
        stimSequence.seq, [0 dur]);
    
ANGLEMAPS = zeros(Uy,Ux,p.nrepeats, p.nstim);
ABSMAPS = zeros(Uy,Ux,p.nrepeats, p.nstim);
for iStim = 1:4
    
    theseEvents = find(stimSequence.seq == iStim);
    
    for iTr = 1:p.nrepeats
        thisEventV = squeeze(periEventV(theseEvents(iTr),:,:));%omit the 1st dimension
        
        
        Tensor =  svdFrameReconstruct(U, thisEventV);
        %should subtract by mimg or else?
        %do detrend?
        
        %Tensor = Tensor - mean(Tensor,3);
        Tensor = Tensor - svdFrameReconstruct(U,mean(mean(avgPeriEventV,1),3)');
        
        [ComplexMaps, AbsMaps, AngleMaps] = FourierMaps(Tensor, Fs, tf_corrected);
        
        if filterFourierMaps
            realMaps_filt = wiener2(real(ComplexMaps), [filtRF filtRF]);
            imagMaps_filt = wiener2(imag(ComplexMaps), [filtRF filtRF]);
            
            ComplexMaps_filt = realMaps_filt + 1i * imagMaps_filt;
            AngleMaps = angle(ComplexMaps_filt);
        end
        %correction of angle produced by stackset.fouriermaps so that it
        %ranges between [-pi pi] (-pi = earliest, pi = latest)
        AngleMaps = -AngleMaps;
        %         AngleMaps(AngleMaps>2*pi) = AngleMaps(AngleMaps>2*pi) - 2*pi;
        %         AngleMaps(AngleMaps<0) = AngleMaps(AngleMaps<0) + 2*pi;
        
        %shift phase
        %use unwrap?
        %AngleMaps(AngleMaps < minPhase) = AngleMaps(AngleMaps < minPhase) + 2*pi;
        %AngleMaps = AngleMaps - 2*pi;
        
        AngleMaps = AngleMaps + delayPhase;
        AngleMaps(AngleMaps > pi) = AngleMaps(AngleMaps > pi) - 2*pi;
        AngleMaps(AngleMaps < -pi) = AngleMaps(AngleMaps < -pi) + 2*pi;
        
        ANGLEMAPS(:, :, iTr, iStim) = AngleMaps;
        ABSMAPS(:,:,iTr, iStim) = AbsMaps;
    end
end

%% combine opposite directions ... needs revision
nRepeat = 4;%p.nrepeats;
randidx = randi(nRepeat,p.nstim,n_boot);
xMap = (ANGLEMAPS(:,:,randidx(1,:),1) - ANGLEMAPS(:,:,randidx(2,:),2));
yMap = (ANGLEMAPS(:,:,randidx(3,:),3) - ANGLEMAPS(:,:,randidx(4,:),4));

xMapm = nanmedian(xMap,3);
yMapm = nanmedian(yMap,3);
xMapm(~mask)=nan;
yMapm(~mask)=nan;

% convert from phase to visual field
% doublePhaseRange = [-2*pi 2*pi*stimDur/(blankDur+stimDur)];%in theory this should be used but
doublePhaseRange = [-pi pi*stimDur/(blankDur+stimDur)];%in practice this must be correct

%range of stimulus position in visual field
%TODO: retrieve this from .x file
xposRange = [-55 55];%deg
yposRange = [-40 40];

xMapm_vf = diff(xposRange)/diff(doublePhaseRange)*(xMapm - doublePhaseRange(1)) + xposRange(1); %[deg]
yMapm_vf = diff(yposRange)/diff(doublePhaseRange)*(yMapm - doublePhaseRange(1)) + yposRange(1); %[deg]

%% compute visual field sign map from xMap and yMap
signMap_boot = nan(size(U,1),size(U,2),n_boot);
for curr_boot = 1:n_boot
    signMap_boot(:,:,curr_boot) = VFS(xMap(:,:,curr_boot), yMap(:,:,curr_boot));
end
signMap_median = imgaussfilt(nanmedian(signMap_boot,3),1);


%% visualize
figname = [dat.constructExpRef(expt.subject, thisDate, thisSeries, expt.expNum) '_'...
    p.xfile(1:end-2) '_' params.movieSuffix];

if params.useCorrected
    figname = [figname '_corrected'];
end

%todo: single trial response at specified pixel as in mainScript_wf.m

figure('position',[0 0 1500 400]);
subplot(131);imagesc(xMapm_vf);axis equal tight; 
% caxis(xposRange);
caxis(prctile(xMapm_vf(:),[1 99]));
mcolorbar(gca,.5);
grid on

subplot(132);imagesc(yMapm_vf);axis equal tight; 
% caxis(yposRange);
caxis(prctile(yMapm_vf(:),[1 99]));
mcolorbar(gca,.5);
grid on

subplot(133);imagesc(signMap_median,'alphadata',mimg./prctile(mimg(:),99));axis equal tight; 
caxis([-1 1]);mcolorbar(gca,.5);
grid on
colormap(jet);

screen2png(fullfile(resultSaveDir,[figname 'signMap.png']));
save(fullfile(resultSaveDir,[figname 'signMap']), ...
    'signMap_median','xMapm','yMapm','xMapm_vf','yMapm_vf','resizeS');
close;




for iStim = 1:p.nstim
    for iRepeat = 1:p.nrepeats
        subplot(p.nstim, p.nrepeats, iRepeat + p.nrepeats*(iStim-1))
        imagesc(squeeze(ANGLEMAPS(:,:,iRepeat,iStim)));
        if iRepeat == 1
            ylabel(['stim ' num2str(iStim)]); 
        end
        if iStim == 1
            title(['repeat ' num2str(iRepeat)]);
        end
    end
end
caxis([-pi pi]);
mcolorbar(gcf, .5);
colormap(hsv);
screen2png(fullfile(resultSaveDir,[figname 'angleMaps.png']));
close;


for iStim = 1:p.nstim
    for iRepeat = 1:p.nrepeats
        subplot(p.nstim, p.nrepeats, iRepeat + p.nrepeats*(iStim-1))
        imagesc(squeeze(ABSMAPS(:,:,iRepeat,iStim)));
        if iRepeat == 1
            ylabel(['stim ' num2str(iStim)]); 
        end
        if iStim == 1
            title(['repeat ' num2str(iRepeat)]);
        end
    end
end
mcolorbar(gcf, .5);
colormap(jet);
screen2png(fullfile(resultSaveDir,[figname 'absMaps.png']));
close;