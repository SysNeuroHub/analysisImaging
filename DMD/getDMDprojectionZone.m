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

% camImg = cameraImageInfo([height width], bregma, lambda, MmPerPixel, binning, topLeftPix);

resized = imresize(DMDfullOn, 1/camImg.binning);

rowIdx = camImg.topLeftPix(1):camImg.topLeftPix(1)+camImg.imageSize(1)-1;
colIdx = camImg.topLeftPix(2):camImg.topLeftPix(2)+camImg.imageSize(2)-1;

cropped = resized(rowIdx, colIdx);

DMDprojectionZone = (cropped>.5*median(cropped(:)));

if showFigure
    figure;
    imagesc(resized);hold on;
    rectangle('Position',[camImg.topLeftPix(2) camImg.topLeftPix(1) camImg.imageSize(2) camImg.imageSize(1)])
    axis equal tight; grid minor;
    title(latestImageFile);
end