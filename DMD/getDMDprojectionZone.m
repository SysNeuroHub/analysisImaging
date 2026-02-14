function DMDprojectionZone = getDMDprojectionZone(camImg, imgDir, latestImageFile, showFigure)
% DMDprojectionZone = getDMDprojectionZone(camImg, imgDir, latestImageFile)

if nargin < 4
    showFigure = 0;
end
if nargin < 2 || isempty(imgDir)
    imgDir = '/mnt/dshi0006_market/DMDprojectionZone';
end
if nargin < 1
    camImg = cameraImageInfo; %'default' position
end

if nargin < 3 || isempty(latestImageFile)
    imageFiles = dir(fullfile(imgDir, '*.tif'));

    % Sort files by date
    [~, idx] = sort([imageFiles.datenum], 'descend');
    latestImageFile = imageFiles(idx(1)).name;
end

DMDfullOn = imread(fullfile(imgDir, latestImageFile));
% DMDfullOn=loadTiffStack('/DMDprojection20260209.tif','tiffobj',1);

cropped = resizeCropFullImg(DMDfullOn, camImg);

DMDprojectionZone = (cropped>.5*median(cropped(:)));

if showFigure
    figure;
    imagesc(imresize(DMDfullOn, 1/camImg.binning));hold on;
    rectangle('Position',[camImg.topLeftPix(2) camImg.topLeftPix(1) camImg.imageSize(2) camImg.imageSize(1)])
    axis equal tight; grid minor;
    title(latestImageFile);
end