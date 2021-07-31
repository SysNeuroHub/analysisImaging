
% tried to obtain RF of each pixel by reverse correlation
%cf. quickPixelMapRet.m by SF and AP_sparse_noise_retinotopy
addpath(genpath('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\master\visbox\Stimulus'));%stimScreen
addpath('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\master\stimFiles\current');%xfiles
addpath(genpath('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\master\Analysis'));
rmpath(genpath('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\master\Analysis\wfAnalysis\obs'));

%% experiment info
expt.subject = '97138';%'dummy_wf';
expt.expDate = '2020-07-03_1';%'2020-06-25_1';
expt.expNum = 1;

%NG SF's stim on this code

%% analysis info
nSV=1000;
params.movieSuffix = 'blue';% 'purple''corr_dFF'; %'blue'  %U for corr_dFF can have NANs..
params.useCorrected = 1;
respWin = [0.1 0.5];%[0.1 0.5];
%AP: 6s = [0.05 0.2], 6f = [0 0.18];
%SF: 6s = [0.1,0.5], 6f = [0.05,0.2]
% blWin = [-1 0]; %SF
resizeS = 0.25;

% For ridge regression
lambda = 10; 

%for VFS
n_boot = 1;%must be one when using RF
use_method = 'max'; % max or com
screen_resize_scale = 3; %3 if max method

resultSaveDir = fullfile('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\wf',...
    [expt.subject '_' expt.expDate '_' num2str(expt.expNum)]);
mkdir(resultSaveDir);

%% data location
thisDate = expt.expDate(1:10);
thisSeries = str2num(expt.expDate(12:end));
expPath = fileparts(dat.expFilePath(expt.subject, thisDate, thisSeries, expt.expNum, 'widefield','master'));
saveVpath = expPath;

%% get stimulus info
expt = grabStimTimesWF(expt,1); %expt.stimTimes
expt.timelineFile = dat.expFilePath(expt.subject, thisDate, thisSeries, expt.expNum, 'timeline','master');
expt = grabSparseFrames(expt);  %expt.stimFrames

%below holds for stimSparseNoiseUncorrAsync (AP) & stimSparseNoiseUncorr2 (SF)
pxSize = expt.stimFrames.ss.Parameters(7)/10; 
aziDeg = expt.stimFrames.ss.Parameters(2)/10:pxSize:expt.stimFrames.ss.Parameters(3)/10;
elvDeg = expt.stimFrames.ss.Parameters(4)/10:pxSize:expt.stimFrames.ss.Parameters(5)/10;


stimInd = find([1 diff(expt.stimFrames.ss.ImageSequence) > 0]);
stimImg = nan(size(expt.stimFrames.ss.ImageTextures{1},1),...
    size(expt.stimFrames.ss.ImageTextures{1},2));
for t = 1:length(stimInd)
    stimImg(:,:,t) = expt.stimFrames.ss.ImageTextures{expt.stimFrames.ss.ImageSequence(stimInd(t))};
end
stimImgFlat = reshape(stimImg,[],size(stimImg,3));   % Change stim to [visPosition(=nyStim*nxStim) stimEvents]

%< upto this point, every experiment in 2p/wf should be able to pass


%% check case for mismatch between photodiode and stimulus in Protocol
frameTimes_pd = expt.stimTimes.frameTimes{1};


if strcmp(expt.stimFrames.ss.Type, 'stimSparseNoiseUncorrAsync') %AP
    stim_screen = cat(3, expt.stimFrames.ss.ImageTextures{:});
elseif strcmp(expt.stimFrames.ss.Type, 'stimSparseNoiseUncorr2') %SF
   stim_screen = nan(size(expt.stimFrames.ss.ImageTextures{1},1),size(expt.stimFrames.ss.ImageTextures{1},2),expt.stimFrames.ss.nFrames);
    for t = 1:expt.stimFrames.ss.nFrames
        stim_screen(:,:,t) = expt.stimFrames.ss.ImageTextures{expt.stimFrames.ss.ImageSequence(t)};
    end
else
    error([expt.stimFrames.ss.Type ' is currently not compatible.']);
end
nyStim = size(stim_screen,1);
nxStim = size(stim_screen,2);

% odd number of stimuli, but one extra photodiode flip to come back down
if mod(size(stim_screen,3),2) == 1 && ...
        length(frameTimes_pd) == size(stim_screen,3) + 1
    frameTimes_pd(end) = [];
    warning('Odd number of stimuli, removed last photodiode');
end

% If there's still a mismatch, break
if size(stim_screen,3) ~= length(frameTimes_pd)
    warning([num2str(size(stim_screen,3)) ' stimuli, ', ...
        num2str(length(frameTimes_pd)) ' photodiode pulses']);
    
    % Try to estimate which stim were missed by time difference
    photodiode_diff = diff(frameTimes_pd);
    max_regular_diff_time = prctile(diff(frameTimes_pd),99);
    skip_cutoff = max_regular_diff_time*1.5;%2;
    photodiode_skip = find(photodiode_diff > skip_cutoff);
    est_n_pulse_skip = ceil(photodiode_diff(photodiode_skip)/max_regular_diff_time)-1;
    stim_skip = cell2mat(arrayfun(@(x) photodiode_skip(x):photodiode_skip(x)+est_n_pulse_skip(x)-1, ...
        1:length(photodiode_skip),'uni',false));
    
    if isempty(est_n_pulse_skip) || length(frameTimes_pd) + sum(est_n_pulse_skip) ~= size(stim_screen,3)
        error('Can''t match photodiode events to stimuli')
    end
end

stim_times = frameTimes_pd;


%% retrieve imaging data (SVD format)
%[U,V,t] = loadSVD(expt);%not yet
[U, V, expt.frameTimes, mimg] = quickLoadUVt(expPath, nSV, saveVpath, params);
Ux0 = size(U,2);
Uy0 = size(U,1);
U = imresize(U,resizeS);
U = U(:,:,1:nSV);
%mimg = imresize(mimg,resizeS);
%movieWithTracesSVD(U,V,expt.frameTimes) 
framerate = 1./nanmedian(diff(expt.frameTimes));

%% some post processing of imaging data (optional)
% subtract global response?

%% stimulus-triggered V
%% AP:
% Get average response to each stimulus
%INPUT: stim_screen, V, respWin
%OUTPUT: response_grid{px_y,px_x}(nSV, ?)
surround_samplerate = 1/(framerate*1);
surround_time = respWin(1):surround_samplerate:respWin(2);
response_n = nan(nyStim,nxStim);
response_grid = cell(nyStim,nxStim);
response_t = [];
%for px_y = 1:nyStim
%    for px_x = 1:nxStim
        
        % Use first frame of dark or light stim
        %align_stims = (stim_screen(:,:,2:end)~= 0) & ...
        %    (diff(stim_screen,[],3) ~= 0);
        %align_times = stim_times(find(align_stims)+1);
        align_times = stim_times(stimInd);
%         keyboard
%         align_times = align_times(round(length(align_times)/2):end);%why do this??
        
%         response_n(px_y,px_x) = length(align_times);
%         keyboard
        % Don't use times that fall outside of imaging
        align_times(align_times + surround_time(1) < expt.frameTimes(2) | ...
            align_times + surround_time(2) > expt.frameTimes(end)) = [];
%         keyboard
        % Get stim-aligned responses, 2 choices:
        
        % 1) Interpolate times (slow - but supersamples so better)
        %         align_surround_times = bsxfun(@plus, align_times, surround_time);
        %         peri_stim_v = permute(mean(interp1(expt.frameTimes,V',align_surround_times),1),[3,2,1]);
        
        % 2) Use closest frames to times (much faster - not different)
        align_surround_times = bsxfun(@plus, align_times, surround_time);
        frame_edges = [expt.frameTimes,expt.frameTimes(end)+1/framerate];
        align_frames = discretize(align_surround_times,frame_edges);
        
        align_frames(any(isnan(align_frames),2),:) = [];
%         keyboard
        % Define the peri-stim V's as subtracting first frame (baseline)
        peri_stim_v = bsxfun(@minus, ...
            reshape(V(:,align_frames)',size(align_frames,1),size(align_frames,2),[]), ...
            reshape(V(:,align_frames(:,1))',size(align_frames(:,1),1),size(align_frames(:,1),2),[]));
%         keyboard
        mean_peri_stim_v = permute(mean(peri_stim_v,2),[3,1,2]); %[nSV, #presentation?] temporal mean of peri_stim_v 
%         keyboard
        % Save V's
%         response_grid{px_y,px_x} = mean_peri_stim_v;
       
        response_t = cat(3,response_t, squeeze(mean(peri_stim_v,1))');
    %end
%end
mV_t = mean(response_t,3);
mresponse_t = squeeze(svdFrameReconstruct(nanmean(nanmean(U)), mV_t));
%somethins is wrong on this mean reaponse

trialResp = mean_peri_stim_v';%correct?
revCor = rReg(stimImgFlat', trialResp, lambda, 0); % Ridge regression [nyStim*nxStim+1, nSV]
revCor3D = reshape(revCor(2:end,:),nyStim,nxStim,nSV);

%pixelMapRFs = spatialQRT(:,iDims)*SVQRT(iDims,iDims)*revCor'; % Get pixel RFs
pixelMapRFs = svdFrameReconstruct(U,revCor'); 
pixelMapRFs = reshape(pixelMapRFs,size(U,1)*size(U,2), nyStim*nxStim+1);%correct??
planeRFs = permute(reshape(pixelMapRFs(:,2:end),[],nyStim,nxStim),[2,3,1]); % [nyStim, nxStim, plane pixel]
%< this already accounted for all stimulus events
planeRFs4D = reshape(planeRFs,size(planeRFs,1),size(planeRFs,2),size(U,1),size(U,2));% [nyStim, nxStim, pixelY pixelX]
%images(squeeze(planeRFs4D(:,:,26:35,30)));
pixelY = 26:35; pixelX=30;
for ii = 1:10
    subplot(1,10,ii);
    imagesc(squeeze(planeRFs4D(:,:,pixelY(ii),pixelX)));
    title(['x:' num2str(pixelX) ', y:' num2str(pixelY(ii))]);
    axis equal tight;
end
caxes(gcf,[0 100],'indirect');mcolorbar;colormap('default');
figname = [dat.constructExpRef(expt.subject, thisDate, thisSeries, expt.expNum) '_'...
    'RF_' params.movieSuffix '_' num2str(params.useCorrected)];
screen2png(fullfile(resultSaveDir,figname),gcf);
close;


figure('position',[0 0 1900 1200]);
subplot(131);plot(surround_time, mresponse_t);xlabel('peri stim time [s]');

%is the below correct??
response_grid = cell(nyStim,nxStim);
for px_y = 1:nyStim
    for px_x = 1:nxStim
        response_grid{px_y, px_x}(:,1) = revCor3D(px_y,px_x,:); 
    end
end

%% compute preferred stimulus in each pixel. cf. sparseRetinotopy by NS
[xMap, yMap] = prefStimSVD(U, response_grid, screen_resize_scale, n_boot, use_method);
xMapm = nanmedian(xMap,3);
yMapm = nanmedian(yMap,3);


%transform stimulus index to visual angle
xMapm = imresize(xMapm, [Uy0 Ux0]);
yMapm = imresize(yMapm, [Uy0 Ux0]);
aziMap = (max(aziDeg)-min(aziDeg))/(length(aziDeg)-1)*(xMapm-1)+aziDeg(1);
elvMap = (max(elvDeg)-min(elvDeg))/(length(elvDeg)-1)*(yMapm-1)+elvDeg(1);

subplot(132);imagesc(aziMap,'alphadata',mimg./prctile(mimg(:),99)); 
set(gca,'color','k');axis equal tight; colorbar;title('pref azimuth');
caxis([min(aziDeg) max(aziDeg)]);caxis([-20 60]);
subplot(133);imagesc(elvMap,'alphadata',mimg./prctile(mimg(:),99)); 
set(gca,'color','k');axis equal tight; colorbar;title('pref elevation');
caxis([min(elvDeg) max(elvDeg)]);caxis([-20 10]);
colormap(jet);

figname = [dat.constructExpRef(expt.subject, thisDate, thisSeries, expt.expNum) '_'...
    'sparseRetinotopy_revCor_test_' params.movieSuffix '_' num2str(params.useCorrected)];
screen2png(fullfile(resultSaveDir,figname),gcf);
