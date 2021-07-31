
%specifications
% needs to be fast so can be used during experiment
% work on SVD data
% analyze response to sparse noise
% compatible with 2p (map and trace) and wf
% options for post processing on imaging data such as dF/F, spikes


addpath(genpath('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\master\visbox\Stimulus'));%stimScreen
addpath('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\master\stimFiles\current\xfiles');%xfiles
addpath(genpath('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\master\Analysis'));
rmpath(genpath('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\master\Analysis\wfAnalysis\obs'));

%% experiment info
expt.subject = 'TIGRE2GCaMP6s_318';
expt.expDate = '2020-12-06_1';
expt.expNum = 2;%9;

%NG SF's stim on this code

%% analysis info
params.modality = 'Fcell';%'FcellCorrected';
params.useCells = 1;
medFiltOrder = 20; %if empty, dont apply median filter
respWin = [0.1 0.5];
%AP: 6s = [0.05 0.2], 6f = [0 0.18];
%SF: 6s = [0.1,0.5], 6f = [0.05,0.2]
%blWin = [-1 0]; %SF

%for estimation of preferred stim
n_boot = 10;%1 to see retinotopy only, 10 to compute VFS
use_method = 'max'; % max or com
screen_resize_scale = 3; %3 if max method

resultSaveDir = fullfile('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\AnalysisResult\2p',...
    [expt.subject '_' expt.expDate '_' num2str(expt.expNum)]);
mkdir(resultSaveDir);

%% data location
thisDate = expt.expDate(1:10);
thisSeries = str2num(expt.expDate(12:end));
% expPath = fileparts(dat.expFilePath(expt.subject, thisDate, thisSeries, expt.expNum, '2p','master'));
% saveS2Ppath = expPath;
marketDir = 'M:\Subjects';

%% get stimulus info
p = ProtocolLoad_wf(expt.subject,expt.expDate,expt.expNum);%10/7/20 now this works if data is uploaded to Market!

expt = grabStimTimes2p(expt,1, marketDir); %expt.stimTimes %SLOW!
%expt.timelineFile = dat.expFilePath(expt.subject, thisDate, thisSeries, expt.expNum, 'timeline','master');
expt = grabSparseFrames(expt);  %expt.stimFrames

%below holds for stimSparseNoiseUncorrAsync (AP) & stimSparseNoiseUncorr2 (SF)
pxSize = expt.stimFrames.ss.Parameters(7)/10; 
aziDeg = expt.stimFrames.ss.Parameters(2)/10:pxSize:expt.stimFrames.ss.Parameters(3)/10;
elvDeg = expt.stimFrames.ss.Parameters(4)/10:pxSize:expt.stimFrames.ss.Parameters(5)/10;


stimInd = find([1 diff(expt.stimFrames.ss.ImageSequence) > 0]);
% stimTimes = stimTimes(stimInd); %SF
stimImg = nan(size(expt.stimFrames.ss.ImageTextures{1},1),...
    size(expt.stimFrames.ss.ImageTextures{1},2));
for t = 1:length(stimInd)
    stimImg(:,:,t) = expt.stimFrames.ss.ImageTextures{expt.stimFrames.ss.ImageSequence(stimInd(t))};
end



%% prepare stim_screen from stimulus protocol
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
ny = size(stim_screen,1);
nx = size(stim_screen,2);

%% check case for mismatch between photodiode and stimulus in Protocol
expt = checkStimFrames(expt, size(stim_screen,3));%10/7/20
%< upto this point, every experiment in 2p/wf should be able to pass



%% retrieve 2p trace ... very slow
[traces, expt.frameTimes, mimg, roi] = quickLoadS2Pt(expt, marketDir, p, params);
nplanes = size(expt.frameTimes,1);
for iplane = 1:nplanes
    
    %movieWithTracesSVD(U,V,expt.frameTimes)
    % framerate = 1./nanmedian(diff(expt.frameTimes));
    
    %% some post processing of imaging data (optional)
    if ~isempty(medFiltOrder)
        traces_cache = cat(2,fliplr(traces{iplane}),traces{iplane}, fliplr(traces{iplane}));
        traces_cache = medfilt1(traces_cache, medFiltOrder,[],2);
        nframes_cache = size(traces{iplane},2);
        traces{iplane} = traces_cache(:,nframes_cache+1:2*nframes_cache);
        clear traces_cache
    end
    
    %% compute retinotopy and visualize it
    ncells = size(traces{iplane},1);
    
    iclust = sum(roi{iplane} .* reshape(1:ncells, 1, 1, ncells), 3); %~dat.res.iclust
    %TODO: deal with cell overwrap
    
    response_grid = getSparseResponse(traces{iplane}, expt.frameTimes(iplane,:), stim_screen, ...
        expt.stimTimes.frameTimes, respWin);
    U = reshape(eye(ncells),ncells,1,ncells); %dummy U to use prefStimSVD and svdFrameReconstruct
    
    %% response map at specified cell in visual field
    figure('position',[0 0 1900 1200]);
    theseCells = 1:ncells;
    nRows = ceil(sqrt(ncells+1));
    nCols = ceil((ncells+1)/nRows);

    response_s = [];
    for icell = 1:length(theseCells)
        for px_y = 1:ny
            for px_x = 1:nx
                response_s(px_y,px_x,icell) = mean(squeeze(svdFrameReconstruct(U(theseCells(icell),1,:),response_grid{px_y,px_x})));%mean across presentations
            end
        end
        subplot(nRows,nCols,icell);
        imagesc(aziDeg, elvDeg, response_s(:,:,icell)); %TODO, translate axes to visual angle
        title(['cell:' num2str(theseCells(icell))]);
        axis equal tight;
    end
    subplot(nRows,nCols,icell+1);
     imagesc(aziDeg, elvDeg,  mean(response_s,3));
     title('avg across cells');
     axis equal tight;
     
     %caxes(gcf,[0 100],'indirect');mcolorbar;colormap('default');
    figname = [dat.constructExpRef(expt.subject, thisDate, thisSeries, expt.expNum) '_'...
        'respMap_' params.modality '_iplane' num2str(iplane)];
    screen2png(fullfile(resultSaveDir,figname),gcf);
    close;
    
    %% FoV modified from roiOnMimg.m
    fig=figure;
    imagesc(mimg(:,:,iplane));colormap(gray);hold on;
    contour(sum(roi{iplane},3),'color',[.9 .4 .4]);
    axis equal tight
    
    for icell = 1:ncells
        c = regionprops(roi{iplane}(:,:,icell),'centroid');
        text(c.Centroid(1),c.Centroid(2),num2str(icell),'color','g');
    end
    
    screen2png(fullfile(resultSaveDir, ['ROI_iplane' num2str(iplane)]));
    
    close(fig);

end

% %% global response across pixels... maynot be correct
% mV_t = mean(response_t,3);
% mresponse_t = squeeze(svdFrameReconstruct(nanmean(nanmean(U)), mV_t));
% figure('position',[0 0 1900 1200]);
% subplot(131);plot(surround_time, mresponse_t);xlabel('peri stim time [s]');

%% compute preferred stimulus in each pixel. cf. sparseRetinotopy by NS
[xMap, yMap] = prefStimSVD(U, response_grid, screen_resize_scale, n_boot, use_method);
xMapm = nanmedian(xMap,3);
yMapm = nanmedian(yMap,3);

%% the code works down to here

% %transform stimulus index to visual angle
% xMapm = nanimresize(xMapm, [Uy0 Ux0]);
% yMapm = nanimresize(yMapm, [Uy0 Ux0]);
% aziMap = (max(aziDeg)-min(aziDeg))/(length(aziDeg)-1)*(xMapm-1)+aziDeg(1);
% elvMap = (max(elvDeg)-min(elvDeg))/(length(elvDeg)-1)*(yMapm-1)+elvDeg(1);
% 
% subplot(132);imagesc(aziMap,'alphadata',mimg./prctile(mimg(:),99)); 
% set(gca,'color','k');axis equal tight; colorbar;title('pref azimuth');
% caxis([min(aziDeg) max(aziDeg)]);caxis([-60 50]);
% subplot(133);imagesc(elvMap,'alphadata',mimg./prctile(mimg(:),99)); 
% set(gca,'color','k');axis equal tight; colorbar;title('pref elevation');
% caxis([min(elvDeg) max(elvDeg)]);caxis([-50 30]);
% colormap(jet);


%% 
%TODO: make 2p version of pixelRFViewer_test
%TODO: dont use dat
iplane = 1;
 s2pname = sprintf('%s/%s/%s/%d/F_%s_%s_plane%d_proc.mat', marketDir, ...
        expt.subject, expt.expDate, expt.expNum, expt.subject, expt.expDate, iplane);
    load(s2pname,'dat');
    
xMapm_c = zeros(1,length(dat.stat));
xMapm_c(find([dat.stat.iscell]))=xMapm;
xMapmn = xMapm_c / nx;
pickCell_test(dat, xMapmn); 
%pickCell_test(dat, [xMapmn yMapmn]); %not yet
%pickCell_test(dat, property, reliability); %not yet 


figname = [dat.constructExpRef(expt.subject, thisDate, thisSeries, expt.expNum) '_'...
    'sparseRetinotopy_' params.modality '_' num2str(params.useCorrected)];
screen2png(fullfile(resultSaveDir,figname),gcf);
