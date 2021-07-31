
%this is a test for SF's stimulus protocol 'stimSparseNoiseUncorr2' for
%pixelRetinotopy with reverse correlation

addpath(genpath('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\master\visbox\Stimulus'));%stimScreen
addpath('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\master\stimFiles\current');%xfiles
addpath(genpath('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\master\Analysis'));
rmpath(genpath('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\master\Analysis\wfAnalysis\obs'));

%% experiment info
expt.subject = '97138';%'dummy_wf';
expt.expDate = '2020-07-03_3';%'2020-06-25_1';
expt.expNum = 2;


%% analysis info
nSV=1000;
params.movieSuffix = 'blue';% 'purple''corr_dFF'; %'blue'  %U for corr_dFF can have NANs..
params.useCorrected = 1;
respWin = [0.1 0.5];%[0.1 0.5];
%AP: 6s = [0.05 0.2], 6f = [0 0.18];
%SF: 6s = [0.1,0.5], 6f = [0.05,0.2]
blWin = [-1 0]; %SF
resizeS = 0.25;

% For ridge regression
lambda = 10; 

% NOTE: Only adjust this if you use a full field stimulus (ALL SCREENS). If
% you adjusted this in mpep, make sure visSpace is [0 1; 0 1]
visSpace = [0 0.75; 0 1]; % Proportion of azimuth, proportion of elevation, 0 = left/top, 1 = right/bottom of stimulus. 

% Set this depending on what side of the brain you are imaging (e.g. left
% visual field for right hemisphere). Including entire stimulus probably
% introduces a higher chance of spurious fits, but you also want to cover the
% entire retinotopic range of the region you imaged

% %for VFS
% n_boot = 10;%1 to see retinotopy only, 10 to compute VFS
% use_method = 'max'; % max or com
% screen_resize_scale = 3; %3 if max method

resultSaveDir = fullfile('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\wf',...
    [expt.subject '_' expt.expDate '_' num2str(expt.expNum)]);
mkdir(resultSaveDir);

%% data location
thisDate = expt.expDate(1:10);
thisSeries = str2num(expt.expDate(12:end));
expPath = fileparts(dat.expFilePath(expt.subject, thisDate, thisSeries, expt.expNum, 'widefield','master'));
saveVpath = expPath;

%% get stimulus info
expt.timelineFile = dat.expFilePath(expt.subject, thisDate, thisSeries, expt.expNum, 'timeline','master');
expt = grabStimTimesWF(expt,1,fileparts(expt.timelineFile),1/5,0); %expt.stimTimes
expt = grabSparseFrames(expt);  %expt.stimFrames

%below holds for stimSparseNoiseUncorrAsync (AP) & stimSparseNoiseUncorr2 (SF)
pxSize = expt.stimFrames.ss.Parameters(7)/10; 
aziDeg = expt.stimFrames.ss.Parameters(2)/10:pxSize:expt.stimFrames.ss.Parameters(3)/10;
elvDeg = expt.stimFrames.ss.Parameters(4)/10:pxSize:expt.stimFrames.ss.Parameters(5)/10;


stimInd = find([1 diff(expt.stimFrames.ss.ImageSequence) > 0]);%stim frame numbers when stimulus changed 
stimImg = nan(size(expt.stimFrames.ss.ImageTextures{1},1),...
    size(expt.stimFrames.ss.ImageTextures{1},2));
for t = 1:length(stimInd)
    stimImg(:,:,t) = expt.stimFrames.ss.ImageTextures{expt.stimFrames.ss.ImageSequence(stimInd(t))};
end



%% check case for mismatch between photodiode and stimulus in Protocol
frameTimes_pd = expt.stimTimes.frameTimes{1};%17995


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

%< upto this point, every noise procol should be able to pass

%% retrieve imaging data (SVD format)
%[U,V,t] = loadSVD(expt);%not yet
[U, V, expt.frameTimes, mimg] = quickLoadUVt(expPath, nSV, saveVpath, params);
Ux0 = size(U,2);
Uy0 = size(U,1);
U = imresize(U,resizeS);
U = U(:,:,1:nSV);
[xRS, yRS, ~] = size(U);
%mimg = imresize(mimg,resizeS);
%movieWithTracesSVD(U,V,expt.frameTimes) 
framerate = 1./nanmedian(diff(expt.frameTimes));


%SF:
stimTimes = frameTimes_pd(stimInd); %TL time when stimulus changes

planeTimes = expt.frameTimes(1:end);

% Remove commonly dropped last frame(s)
if length(planeTimes) > size(V,2)
    lenDiff = length(planeTimes) - size(V,2);
    planeTimes(end-lenDiff+1:end) = [];
end

respTimes = [stimTimes + respWin(1) stimTimes + respWin(2)];
baseTimes = [stimTimes + blWin(1) stimTimes + blWin(2)];

 % Get trial responses for temporal component

 % Apparentally these lines don't work in older versions of matlab... at
 % least 2016a
 try
     respInd = planeTimes >= respTimes(:,1) & planeTimes <= respTimes(:,2);
     baseInd = planeTimes >= baseTimes(:,1) & planeTimes <= baseTimes(:,2);
 catch % These lines should work in older versions
     respInd = bsxfun(@ge, planeTimes, respTimes(:,1)) & bsxfun(@le, planeTimes, respTimes(:,2));
     baseInd = bsxfun(@ge, planeTimes, baseTimes(:,1)) & bsxfun(@le, planeTimes, baseTimes(:,2));
 end

resp = respInd*V'; % Sum activity in resp window for each trial
base = baseInd*V'; % Sum activity in baseline window for each trial
try % Same code issue as above
    trialResp = resp./sum(respInd,2)-base./sum(baseInd,2); % subtract mean of baseline windows from mean of resp windows
catch
    trialResp = bsxfun(@rdivide, resp, sum(respInd,2)) - bsxfun(@rdivide, base, sum(baseInd,2));
end
%trialResp: [all stimulus events, nSV]

stimImgFlat = reshape(stimImg,[],size(stimImg,3));   % Change stim to [visPosition(=nyStim*nxStim) stimEvents]

revCor = rReg(stimImgFlat', trialResp, lambda, 0); % Ridge regression of V [nyStim*nxStim+1, nSV]
% revCor3D = reshape(revCor(2:end,:),nyStim,nxStim,nSV);

%pixelMapRFs = spatialQRT(:,iDims)*SVQRT(iDims,iDims)*revCor'; % Get pixel RFs
pixelMapRFs = svdFrameReconstruct(U,revCor'); 
pixelMapRFs = reshape(pixelMapRFs,size(U,1)*size(U,2), nyStim*nxStim+1);%correct??
planeRFs = permute(reshape(pixelMapRFs(:,2:end),[],nyStim,nxStim),[2,3,1]); % [nyStim, nxStim, plane pixel]
%< this already accounted for all stimulus events


%% Find max of pixel RF maps for retinotopy 

smWidth = 2;

planeRet = size(planeRFs,3);

tic
for p = 1:size(planeRFs,3)
    
    tmp = imgaussfilt(planeRFs(:,:,p),smWidth); 
    [~,mInd] = max(tmp(:));                              
    [x, y] = ind2sub(size(tmp),mInd);
    
    planeRet(p) = x + y*1i; % Save x and y coords as complex number
    
end

planeRet = reshape(planeRet, xRS, yRS); % Reshape to 2d imag

toc
disp('Pixel retinotopy complete.')



%% Plot plane RF maps for azimuth and elevation

% pxSize = round((270*diff(visSpace(1,:)))/nxStim,1);
pxSize = expt.stimFrames.ss.Parameters(7)/10;

aziDeg = expt.stimFrames.ss.Parameters(2)/10:pxSize:expt.stimFrames.ss.Parameters(3)/10;
elvDeg = expt.stimFrames.ss.Parameters(4)/10:pxSize:expt.stimFrames.ss.Parameters(5)/10;

% Set these if you displayed stim on all screens
aziDeg = aziDeg(round(length(aziDeg)*visSpace(1,1))+1:round(length(aziDeg)*visSpace(1,2))); 
elvDeg = elvDeg(round(length(elvDeg)*visSpace(2,1))+1:round(length(elvDeg)*visSpace(2,2)));
retSm = 2; % Smoothing of pixel map

% Find good color range for retinotopy

azi = imag(planeRet(:));
elv = real(planeRet(:));

muAzi = mean(azi);
stdAzi = std(azi);

muElv = mean(elv);
stdElv = std(elv);

cMinElv = floor(muElv - 2*stdElv);
cMaxElv = ceil(muElv + 2*stdElv);

if cMinElv < 1
    cMinElv = 1;
end

if cMaxElv > nyStim
    cMaxElv = nyStim;
end

cMinAzi = floor(muAzi - 2*stdAzi);
cMaxAzi = ceil(muAzi + 2*stdAzi);

if cMinAzi < 1
    cMinAzi = 1;
end

if cMaxAzi > nxStim
    cMaxAzi = nxStim;
end

% If above color range is off, set it manually

figure('position',[0 0 2400 1500]);

% Plot elevation

subplot(2,2,1)
if retSm > 0
    imagesc(imgaussfilt(real(planeRet), retSm),[cMinElv cMaxElv]);
else
    imagesc(real(planeRet), [cMinElv cMaxElv]);
end

box off
colormap jet
c1 = colorbar;
c1.Ticks = cMinElv:2:cMaxElv;
c1.TickLabels = round(elvDeg(c1.Ticks));
c1.TickLabelsMode = 'manual';
c1.Label.String = 'Visual degrees';
axis equal off
title('Elevation')

% Plot azimuth

subplot(2,2,2)
if retSm > 0
    imagesc(imgaussfilt(imag(planeRet), retSm), [cMinAzi cMaxAzi]);
else
    imagesc(imag(planeRet), [cMinAzi cMaxAzi]);
end

box off
colormap jet
c2 = colorbar;
c2.Ticks = cMinAzi:5:cMaxAzi;
c2.TickLabels = round(aziDeg(c2.Ticks));
c2.TickLabelsMode = 'manual';
c2.Label.String = 'Visual degrees';
axis equal off
title('Azimuth')

% mImg = subplot(2,2,3.5);
% imagesc(mean(reshape(IMG, xRS, yRS, []),3));
% box off
% axis equal off
% title('Mean image')
% colormap(mImg,gray)

screen2png('pixel-base retinotopy');

% % %is the below correct??
% % response_grid = cell(nyStim,nxStim);
% % for px_y = 1:nyStim
% %     for px_x = 1:nxStim
% %         response_grid{px_y, px_x}(:,1) = abs(revCor3D(px_y,px_x,:)); %take abs as this is RF not response
% %     end
% % end
% % 
% % %% compute preferred stimulus in each pixel. cf. sparseRetinotopy by NS
% % [xMap, yMap] = prefStimSVD(U, response_grid, screen_resize_scale, n_boot, use_method);
% % xMapm = nanmedian(xMap,3);
% % yMapm = nanmedian(yMap,3);
% % 
% % % %% visual field sign map
% % % signMap_boot = nan(size(U,1),size(U,2),n_boot);
% % % for curr_boot = 1:n_boot
% % %     signMap_boot(:,:,curr_boot) = VFS(xMap(:,:,curr_boot), yMap(:,:,curr_boot));
% % % end
% % % signMap_median = imgaussfilt(nanmedian(signMap_boot,3),2);
% % 
% % 
% % %transform stimulus index to visual angle
% % xMapm = imresize(xMapm, [Uy0 Ux0]);
% % yMapm = imresize(yMapm, [Uy0 Ux0]);
% % aziMap = (max(aziDeg)-min(aziDeg))/(length(aziDeg)-1)*(xMapm-1)+aziDeg(1);
% % elvMap = (max(elvDeg)-min(elvDeg))/(length(elvDeg)-1)*(yMapm-1)+elvDeg(1);
% % 
% % subplot(132);imagesc(aziMap,'alphadata',mimg./prctile(mimg(:),99)); 
% % set(gca,'color','k');axis equal tight; colorbar;title('pref azimuth');
% % caxis([min(aziDeg) max(aziDeg)]);caxis([-20 60]);
% % subplot(133);imagesc(elvMap,'alphadata',mimg./prctile(mimg(:),99)); 
% % set(gca,'color','k');axis equal tight; colorbar;title('pref elevation');
% % caxis([min(elvDeg) max(elvDeg)]);caxis([-20 10]);
% % colormap(jet);
% % 
% % figname = [dat.constructExpRef(expt.subject, thisDate, thisSeries, expt.expNum) '_'...
% %     'sparseRetinotopy_test_' params.movieSuffix '_' num2str(params.useCorrected)];
% % screen2png(fullfile(resultSaveDir,figname),gcf);
