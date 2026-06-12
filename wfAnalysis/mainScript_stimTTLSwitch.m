%% this is a script to do the following:
% load data from Market server
% detect onset of each stimulation in thorsync time using photodiode signal
% detect wf frame time in thorsync time
% obtain stimulus triggered time courses
% compute Direction selectivity index of each cell



setPath_analysisImaging;

%% experiment

expt.subject = 'Confucious';
expt.expDate = '2026-05-21_1';
expt.expNum = 2;
bklightCtrl = 0;
polyScanPro = false;

params.movieSuffix = 'amber';
params.useCorrected = 0;

%% SVD
nSV = 2000;

%% analysis
marginT = 0.2; %[s]
resizeS = 0.5; %spatial rescaling factor



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
[newU, newV] = dffFromSVD(U, V, mimg); %NG

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
laserTh = 1; %[V]
expt = grabLaserTimesWF(expt,[],laserTh);

%% DMD pattern index at the time of laser onset
nPtn = p.pars(2);
[DMDIn_state, DMDOut_state] = getDMDState(Timeline, nPtn);
%% sanity check DMD pattern idx
[f_sanitycheck, NG_min] = DMD_ptn_sanityCheck(expt, Timeline, stimSequence, DMDIn_state, DMDOut_state);
screen2png([figname '_DMD_ptn_sanityCheck'], f_sanitycheck);
    

%check number of stimulation
if length(expt.stimTimes.onset) ~= p.nstim * p.nrepeats
    error(['Detected #stim onset:' num2str(expt.stimTimes.onset') '/ total #stim' num2str(p.nstim * p.nrepeats)]);
end

%% stimulus triggered movie
% pixelTuningCurveViewerSVD(U, V, t, expt.laserTimes.onset, stimSequence.seq, respWin, 1);
% title(tname);

%% save optogen triggered stim and response 
% load('M:\DMD images\CCFBL_584x450pix_5x8circle_l\CCFBL_584x450pix_5x8circle_l_Confucious.mat'); %TMP
load('M:\DMD images\CCFBL_584x450pix_1x1circle_r\CCFBL_584x450pix_1x1circle_r_Confucious.mat',...
    'image4OI','image4OI_all_wCCF')
image4OI = cat(3, image4OI, zeros(size(image4OI,1),size(image4OI,2))); %HACK to add blank stimuli

if ~polyScanPro
    image4OI = round(image4OI);
end
respWin = [-0.5 1.5];%[-0.5 2.5];
[DMD4OIsequence, tevent_DMD] = getDMD4OI_eventLockedAvg(image4OI, ...
    DMDIn_state, Timeline, t, stimSequence, expt.stimTimes, p.nstim, respWin); %NG

[avgPeriEventV, winSamps, periEventV, sortLabels] = ...
    eventLockedAvg(V, t, expt.stimTimes.onset, stimSequence.seq, respWin);
avgPeriEventV = avgPeriEventV - mean(avgPeriEventV(:,:,winSamps<0),3);

% clim_in = [0 4.5];
% clim_out = [-20 20];
% k2cmap = customcolormap([0 1],{'#00FFFF','#000000'});
% k2amap = customcolormap([0 1],{'#FFA500','#000000'});
% figure('position',[50 50 1800 900]);
% ha = tight_subplot(2,p.nstim,[.01 .01],[.01 .01],[.01 .05]);
% position = [1:12];
% for istim = 1:p.nstim
%     inout{istim} = squeeze(DMD4OIsequence(:,:,:,istim));
%     ax(istim) = ha(position(istim));%subplot(8,5,position(istim));
%     clim{istim} = clim_in;
%     cmap{istim} = k2cmap;
%     alphadata{istim} = image4OI_all_wCCF;
% end
% for istim = 1:p.nstim
%     inout{6+istim} = svdFrameReconstruct(U, squeeze(avgPeriEventV(istim,:,:)));
%     ax(6+istim) = ha(position(6+istim));%subplot(8,5,position(istim));
%     clim{6+istim} = clim_out;
%     cmap{6+istim} = k2amap;
%     alphadata{6+istim} = imresize(1-image4OI_all_wCCF, [size(U,1),size(U,2)]);
% end
% playMatrix(inout,clim,...
%     ax,'colormap',cmap,...
%     'alphadata',alphadata,...
%     'framerate',5,...
%     'timevec',winSamps,'saveVideo',true,...
%     'mvName',[figname '_' num2str(istim)],...
%     'quality',90,'bgBlack',false);


clim_out = [-20 20];
k2amap = customcolormap([0 1],{'#FFA500','#000000'});
figure('position',[50 50 700 900]);
ha = tight_subplot(8,5,[.01 .01],[.05 .05],[.01 .01]);
delete(ha([2 3 6 11 40]));
position = [36:39 31:35 26:30 21:25 16:20 12:15 7:10 4:5 1];
for istim = 1:p.nstim
    inout{istim} = svdFrameReconstruct(U, squeeze(avgPeriEventV(istim,:,:)));
    ax(istim) = ha(position(istim));%subplot(8,5,position(istim));
end
inout = circshift(inout,1);
ax = circshift(ax,1);
cbar = zeros(p.nstim,1);cbar(1)=1;
showtime = cbar;
playMatrix(inout,clim_out,...
    ax,'colormap',k2amap,...
    'framerate',5,...
    'alphadata', imresize(1-image4OI_all_wCCF, [size(U,1) size(U,2)]),...
    'timevec',winSamps,'saveVideo',true,...
    'mvName',[figname '_' num2str(istim)],...
    'quality',90,'bgBlack',false,...
    'colorbar',cbar,'showtime',showtime);

% 
% %% time-avg response & stimulus preference map
% for ii = 1:5
% [avgPeriEventV, winSamps, periEventV, sortLabels] = ...
%     eventLockedAvg(V, t, expt.stimTimes.onset(1:ii*10*6), stimSequence.seq(1:ii*10*6), respWin);
% 
% preIdx = find(winSamps<0);
% %[~,postIdx] = min(abs(winSamps - 0.65));%to examine laser artifact
% postIdx = intersect(find(winSamps>0), find(winSamps < 1));%min(p.pfiledurs)));
% 
% %subtract by prestimulus in each condition
% tavgRespV = mean(avgPeriEventV(:,:,postIdx),3) - mean(avgPeriEventV(:,:,preIdx),3);
% %condition x nSV
% 
% tavgResp = svdFrameReconstruct(U, tavgRespV');
% %tavgResp = tavgResp - tavgResp(:,:,p.blankstims);%still blood vessel remains...
% 
% 
% 
% %% single trial traces
% % figure;
% %     yy=100; xx = 215;
% %     icond = 6;
% %     these = find(sortLabels == icond);
% %     trace = [];
% %     for ii = 1:numel(these)
% %         trace(:,ii)=svdFrameReconstruct(U(yy,xx,:),squeeze(periEventV(these(ii),:,:)));
% %     end
% %     subplot(121);
% %     imagesc(winSamps, 1:numel(these), trace');vline(0);
% %     title('F')
% %     caxis([-27*3 34*3])
% %     mcolorbar(gca, .5);
% %     xlabel('Time from laser onset [s]');
% %     ylabel('trial');
% % 
% %     subplot(122);
% %     imagesc(winSamps, 1:numel(these), trace'-mean(trace(preIdx,:))');vline(0);
% %     title('F-F0')
% %     caxis([-27*3 34*3])
% %     mcolorbar(gca, .5);
% 
%     
% %% avg Resp in dF/F    
% range_c = 100*squeeze(tavgResp./mimg);
% % crange = max(abs(prctile(range_c(:),[1 99.99])));
% crange = [1];
% 
% showImage = 1;
% nRows            = 1;%ceil(sqrt(p.nstim));
% nCols = p.nstim;%ceil(p.nstim/nRows);
% figure('position',[0 0 1900 900]);
% panel = [];
% for istim = 1:p.nstim
%     panel(istim) = subplot(nRows,nCols,istim);
%     if showImage
%         imagesc(100*squeeze(tavgResp(:,:,istim)./mimg),'alphadata',mask);
%         axis equal tight off;
%     else
%         avgResp = svdFrameReconstruct(U,  squeeze(avgPeriEventV(istim,:,:)));
%         savgResp = squeeze(mean(mean(avgResp(171:190, 61:80,:)./mimg(171:190, 61:80),1),2));
%         plot(winSamps, 100*squeeze(savgResp))
%     end
%     ptnIdx = p.pars(3,istim);
%     mV = p.pars(4,istim);
%     stimDurms = p.pars(6,istim) - p.pars(5,istim);
%     title(sprintf('ptn %d\ninput %d[mV]\nstimDur %d[ms]',ptnIdx,mV,stimDurms));
%     
%     if showImage
%         caxis([-crange crange]);
%     else
%         ylim([-crange crange]);
%         vbox(0, 1e-3*stimDurms);
%     end
% end
% [h,g]=mcolorbar(gca,.5);
% %linkaxes(panel);
% %     xlim([40 110]); ylim([140 220])
% g.YLabel.String='dF/F [%]';
% screen2png(fullfile(resultSaveDir,[figname '_' num2str(10*ii) 'repeats']));
% close;
% end
% 


