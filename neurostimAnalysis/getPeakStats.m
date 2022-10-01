function [xp,yp,tp,pTrace,pImage] = getPeakStats(singlePeriEventStack, TimeVec, ...
    stimInfo, fixedLatency, nanMask, pth, colorName)
doMedian = false;
clusterPermutation = false;
bonferroni = false;

if nargin < 6
    pth = 0.01;
end
if nargin < 7
    colorName = @parula;
end

[imageSize_r(1),imageSize_r(2),nFrames, nTrials] = size(singlePeriEventStack);
nConds = numel(unique(stimInfo.stimLabels));

avgPeriEventStack = [];
for icond = 1:nConds
    theseEvents = find(stimInfo.stimLabels == stimInfo.condLabels(icond));
    avgPeriEventStack(:,:,:,icond) = squeeze(mean(singlePeriEventStack(:,:,:,theseEvents),4));
end

yp = zeros(nConds,1); xp=zeros(nConds,1); tp=zeros(nConds,1);
for icond = 1:nConds
    if isempty(fixedLatency)
        Values_tmp = avgPeriEventStack(:,:,:,icond) .* nanMask;
        [yp(icond), xp(icond), tp(icond)] = getTensorPeak(Values_tmp, 'min');
    else
        [~,tp(icond)] = min(abs(TimeVec - fixedLatency));
        Values_tmp = avgPeriEventStack(:,:,tp(icond),icond) .* nanMask;
        [yp(icond), xp(icond)] = getTensorPeak(Values_tmp, 'min');
    end
end
cp = colorName(nConds);


%% single trial analysis in time
meanTrace = zeros(numel(TimeVec),nConds);
seTrace = zeros(numel(TimeVec), nConds);
pTrace = zeros(numel(TimeVec), nConds);
for icond = 1:nConds
    theseEvents = find(stimInfo.stimLabels == stimInfo.condLabels(icond));
    theseTraces = squeeze(singlePeriEventStack(yp(icond),xp(icond),:,theseEvents));
    if doMedian
        meanTrace(:,icond) = squeeze(nanmedian(theseTraces,2));
        seTrace(:,icond) = 1/sqrt(numel(theseEvents)) * squeeze(mad(theseTraces,[],2));%incorrect
    else
        meanTrace(:,icond) = squeeze(nanmean(theseTraces,2));
        seTrace(:,icond) = 1/sqrt(numel(theseEvents)) * squeeze(std(theseTraces,[],2));
    end
    
    if clusterPermutation
        [clusters, p_values, t_sums, permutation_distribution ] = permutest(...
            zeros(size(theseTraces)),theseTraces, false, pth,[],false); %only detect negative deflection
        pTrace(:,icond) = 1;
        pTrace(clusters{1},icond) = 0;
    else
        cellTraces = mat2cell(theseTraces, ones(size(theseTraces,1),1));
        pTrace(:,icond) = cellfun(@signrank, cellTraces);
    end
end
if ~clusterPermutation && bonferroni
    pTrace = numel(pTrace(:))*pTrace;
end


%% single trial analysis in space
pImage = zeros(imageSize_r(1), imageSize_r(2), nConds);
for icond = 1:nConds
    theseEvents = find(stimInfo.stimLabels == stimInfo.condLabels(icond));
    theseImages = squeeze(singlePeriEventStack(:,:,tp(icond),theseEvents));
    
    if clusterPermutation %SLOW!
        [clusters, p_values, t_sums, permutation_distribution ] = permutest(...
            zeros(size(theseImages)),theseImages, false, pth,[],false);
        %only detect negative deflection, 1 cluster
        pImage_c = ones(prod(imageSize_r),1);
        pImage_c(clusters{1}) = 0;
        pImage(:,:,icond) = reshape(pImage_c, imageSize_r(1), imageSize_r(2));
    else
        cellImages = num2cell(theseImages,[3]);
        pImage(:,:,icond) = cellfun(@(x)signrank(squeeze(x)),cellImages);
    end
end
if ~clusterPermutation && bonferroni
    pImage = numel(pImage(:))*pImage;
end


%% images at the peak in each condition
clim = prctile(avgPeriEventStack(:),[0 99]);

figure('position',[0 0 1920 1080]);
for icond = 1:nConds
    subplot(3,nConds,icond)
    imagesc(avgPeriEventStack(:,:,tp(icond),icond));colormap(gray);
    axis equal tight; hold on;
    vline(xp(icond),[],[],cp(icond,:));
    hline(yp(icond),[],[],cp(icond,:));
    if icond==1
    title(sprintf('%s\n%d, delay %.1f', stimInfo.labelDescription, stimInfo.condLabels(icond),TimeVec(tp(icond))));
    else
    title(sprintf('%d, delay %.1f',stimInfo.condLabels(icond),TimeVec(tp(icond))));
    end
    caxis(clim);
    
    subplot(3,nConds,nConds+icond)
    %plot(Savg.TimeVec, squeeze(Savg.Values(yp(icond),xp(icond),:,icond)),'color',cp(icond,:));
    shadedErrorBar(TimeVec, meanTrace(:,icond),seTrace(:,icond));
    [ss,ee] = trace2box(pTrace(:,icond),pth);
    ylim(clim);
    vline(TimeVec(tp(icond)));
    vbox(TimeVec(ss), TimeVec(ee),[],[cp(icond,:) .5])
 
    title([num2str(stimInfo.condLabels(icond)) ', x=' num2str(xp(icond)) ', y=' num2str(yp(icond))]);
end

% subplot(313);
% % imagesc(imageData.meanImage);
% imagesc(zeros(imageSize_r));
% colormap(gray);
% hold on
% for icond = 1:nConds
%     imSuper(gca,pImage(:,:,icond)<pth,cp(icond,:),0.3);
% end
% scatter(xp,yp,20,cp,'filled','markeredgecolor','w');
% axis equal tight;
