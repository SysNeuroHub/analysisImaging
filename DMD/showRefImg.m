function f = showRefImg(camImg, imgDir, latestImageFile, showFigure)
% f = png4DMD_ref(camImg, imgDir, latestImageFile, showFigure)
% creates a reference image with CCF, bregma, lambda, scalebar and DMD
% projection zone. This should be uniquely defined by cameraImageInfo
% to save the resulting figure, use exportPng4DMD(fullfile(saveDir, [saveName_s '_ref']), f, binarise);

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

f=figure;
DMDprojectionZone = getDMDprojectionZone(camImg, imgDir, latestImageFile, showFigure);
image(zeros(camImg.imageSize));
addAllenCtxOutlines(camImg.bregmapix, camImg.lambdapix, 'w', camImg.MmPerPixel);%this looks at lambda and shrinks the CCF
scatter(camImg.bregmapix(2), camImg.bregmapix(1), 600/camImg.binning,'markerEdgecolor','w','MarkerFaceColor','w','marker','x');
scatter(camImg.lambdapix(2), camImg.lambdapix(1), 600/camImg.binning,'markerEdgecolor','w','MarkerFaceColor','w','marker','x');
line([camImg.lambdapix(2) camImg.lambdapix(2)], [1 camImg.imageSize(1)], 'color','w');
line([camImg.imageSize(2)/2-50/camImg.binning camImg.imageSize(2)/2-50/camImg.binning+1/camImg.MmPerPixel], [800 800]/camImg.binning,'linewidth',2,'color','w');
contour(DMDprojectionZone,'w');