
%specifications
% needs to be fast so can be used during experiment
% work on SVD data
% analyze response to sparse noise
% compatible with 2p (map and trace) and wf
% options for post processing on imaging data such as dF/F, spikes

%% TODO
%subtract global signal across pixels?

setPath_analysisImaging;

%% experiment info
expt.subject = 'hercules';
expt.expDate = '2025-05-29_1';%'2020-06-25_1';%'2020-06-16_1';
expt.expNum = 5;%9;
bklightCtrl = 0; %6/11/20

pixies = [50 50];%selected pixels to show RF [

%NG SF's stim on this code

%% analysis info
nSV=1000;
params.movieSuffix = 'red';% 'purple''corr_dFF'; %'blue'  %U for corr_dFF can have NANs..
params.useCorrected = 0;
respWin = [0.1 0.5];
%AP: 6s = [0.05 0.2], 6f = [0 0.18];
%SF: 6s = [0.1,0.5], 6f = [0.05,0.2]
%blWin = [-1 0]; %SF
resizeS = .25; %.125 Resize factor in cortical space

%for estimation of preferred stim
n_boot = 10;%1 to see retinotopy only, 10 to compute VFS
use_method = 'max'; % max or com
screen_resize_scale = 3; %3 if max method

resultSaveDir = fullfile('M:\Subjects\',expt.subject,...
    expt.expDate, num2str(expt.expNum));

%% data location
thisDate = expt.expDate(1:10);
thisSeries = str2num(expt.expDate(12:end));
expPath = fileparts(dat.expFilePath(expt.subject, thisDate, thisSeries, expt.expNum, 'widefield','master'));
saveVpath = expPath;

%% get each stimulation frame
p = ProtocolLoad_wf(expt.subject,expt.expDate,expt.expNum); %3/6/20
expt = grabStimTimesWF(expt,[p.nstim p.nrepeats],[],[],bklightCtrl); %expt.stimTimes
%expt.timelineFile = dat.expFilePath(expt.subject, thisDate, thisSeries, expt.expNum, 'timeline','master');
expt = grabSparseFrames(expt);  %expt.stimFrames

%below holds for stimSparseNoiseUncorrAsync (AP) & stimSparseNoiseUncorr2 (SF)
iStim = 1;
pxSize = expt.stimFrames(iStim).ss.Parameters(7)/10; 
aziDeg = expt.stimFrames(iStim).ss.Parameters(2)/10:pxSize:expt.stimFrames(iStim).ss.Parameters(3)/10;
elvDeg = expt.stimFrames(iStim).ss.Parameters(4)/10:pxSize:expt.stimFrames(iStim).ss.Parameters(5)/10;


% % stimInd = find([1 diff(expt.stimFrames.ss.ImageSequence) > 0]);
% % % stimTimes = stimTimes(stimInd); %SF
% % stimImg = nan(size(expt.stimFrames.ss.ImageTextures{1},1),...
% %     size(expt.stimFrames.ss.ImageTextures{1},2));
% % for t = 1:length(stimInd)
% %     stimImg(:,:,t) = expt.stimFrames.ss.ImageTextures{expt.stimFrames.ss.ImageSequence(stimInd(t))};
% % end
%< upto this point, every experiment in 2p/wf should be able to pass

[stim_screen, ny, nx] = getStimScreen(expt.stimFrames);

% % %% check case for mismatch between photodiode and stimulus in Protocol
nFramesInProtocol=[];
for iStim = 1:p.nstim
    nFramesInProtocol(iStim) = size(stim_screen,3);
end
expt = checkStimFrames(expt, nFramesInProtocol);%10/7/20
%expt.stimTimes.frameTimes = expt.stimTimes.frameTimes([1:3]); %TEMP

%% retrieve imaging data (SVD format)
[U, V, expt.frameTimes, mimg, mask] = quickLoadUVt(expPath, nSV, saveVpath, params);
Ux0 = size(U,2);
Uy0 = size(U,1);
U = imresize(U,resizeS);
U = U(:,:,1:nSV);
mask = imresize(mask,resizeS);
mimg_rs = imresize(mimg,resizeS);
%movieWithTracesSVD(U,V,expt.frameTimes) 
% framerate = 1./nanmedian(diff(expt.frameTimes));

%% some post processing of imaging data (optional)
% subtract global response?

%% compute response to each stimulus grid 10/11/20
[response_grid, surround_time, response_t] = getResponseGridV(expt, stim_screen, V, respWin);
 
%% global response across pixels... maynot be correct
mV_t = mean(response_t,3);
mresponse_t = squeeze(svdFrameReconstruct(nanmean(nanmean(U)), mV_t));
figure('position',[0 0 1900 1200]);
subplot(131);plot(surround_time, mresponse_t);xlabel('peri stim time [s]');

%% compute preferred stimulus in each pixel. cf. sparseRetinotopy by NS
[xMap, yMap] = prefStimSVD(U, response_grid, screen_resize_scale, n_boot, use_method);
xMapm = nanmedian(xMap,3);
yMapm = nanmedian(yMap,3);
xMapm(~mask)=nan;
yMapm(~mask)=nan;

%% visual field sign map
signMap_boot = nan(size(U,1),size(U,2),n_boot);
for curr_boot = 1:n_boot
    signMap_boot(:,:,curr_boot) = VFS(xMap(:,:,curr_boot), yMap(:,:,curr_boot));
end
signMap_median = imgaussfilt(nanmedian(signMap_boot,3),2);


%transform stimulus index to visual angle
xMapm = nanimresize(xMapm, [Uy0 Ux0]);
yMapm = nanimresize(yMapm, [Uy0 Ux0]);
aziMap = (max(aziDeg)-min(aziDeg))/(length(aziDeg)-1)*(xMapm-1)+aziDeg(1);
elvMap = (max(elvDeg)-min(elvDeg))/(length(elvDeg)-1)*(yMapm-1)+elvDeg(1);

subplot(132);imagesc(aziMap,'alphadata',mimg./prctile(mimg(:),99)); 
set(gca,'color','k');axis equal tight; colorbar;title('pref azimuth');
caxis([min(aziDeg) max(aziDeg)]);%caxis([-50 60]);
subplot(133);imagesc(elvMap,'alphadata',mimg./prctile(mimg(:),99)); 
set(gca,'color','k');axis equal tight; colorbar;title('pref elevation');
caxis([min(elvDeg) max(elvDeg)]);%caxis([-50 10]);
colormap(jet);

figname = [dat.constructExpRef(expt.subject, thisDate, thisSeries, expt.expNum) '_'...
    'sparseRetinotopy_test_' params.movieSuffix '_' num2str(params.useCorrected) ...
    '_' num2str(1e2*respWin(1)) '-' num2str(1e2*respWin(2)) 'ms'];
screen2png(fullfile(resultSaveDir,figname),gcf);


%% response map at specified pixel
%TODO: use bootstrap as above to estimate RF
RF = zeros(size(U,1),size(U,2),ny,nx);
for px_y = 1:ny
    for px_x = 1:nx
        RF(:,:,px_y,px_x) = mean(svdFrameReconstruct(U,response_grid{px_y,px_x}),3);%mean across presentations
    end
end

if ~isempty(pixies)
    for ff = 1:size(pixies,1)
        pixelRFViewer_test(mimg_rs, RF, aziDeg, elvDeg, pixies(ff,:));
        screen2png(fullfile(resultSaveDir,[figname ' at' num2str(pixies(ff,:))]),gcf);
        close(gcf);
    end
end

% caxes(gcf,[0 100],'indirect');mcolorbar;colormap('default');
% figname = [dat.constructExpRef(expt.subject, thisDate, thisSeries, expt.expNum) '_'...
%     'respMap_' params.movieSuffix '_' num2str(params.useCorrected)];
% screen2png(fullfile(resultSaveDir,figname),gcf);
% close;

% % %% showing RF of individual cells
% % % TODO: automatic cell pos detection
% % cellPos = [129 195; 149 290; 155 410; 342 90; 259 555; 406 312; 310 443; 360 385];
% % ncells = size(cellPos,1);
% % 
% % for icell = 1:ncells
% %     subplot(2,4,icell);
% %     imagesc(aziDeg, elvDeg, squeeze(RF(cellPos(icell,1),cellPos(icell,2),:,:)));
% %     title(['ROI:' num2str(icell) ', pixel (' num2str(cellPos(icell,1)) ',' num2str(cellPos(icell,2)),')']);
% % end
% % figname = [dat.constructExpRef(expt.subject, thisDate, thisSeries, expt.expNum) '_'...
% %     'sparseRetinotopy_RFs_' params.movieSuffix '_' num2str(params.useCorrected)];
% % screen2png(fullfile(resultSaveDir,figname),gcf);
