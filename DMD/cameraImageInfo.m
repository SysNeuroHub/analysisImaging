function info = cameraImageInfo(imageSize, bregmapix, lambdapix, MmPerPixel, binning, topLeftPix)
%info = cameraImageInfo(imageSize, bregmapix, lambdapix, MmPerPixel, binning, topLeftPix)
% INPUT
%binning: camera binning factor [1, 2, 4] (can be used as a scale factor)
%imageSize: size of the image recorded by the camera, after binning [height, width] in pix
%bregmaPix: location of the breagma [y x] in the image in pix
%lambdaPix: location of the lambda [y x] in the image in pix
%MmPerPixel: millimer/pixel after binning
%topLeftPix: location of the top left of ROI after binning [y x]
%
% NOTE on 7/7/26:
% These information is critical in generating DMD images in stereo coords,
% such as showAlllenCCFBregmaLambda.m and showRefImg
% These information is not critical in aligning widefield image to stereo
% or CCF coords

if nargin < 1 || isempty(imageSize)
    imageSize = [900 1168];
end
if nargin < 2 || isempty(bregmapix)
    bregmapix = [360 imageSize(2)/2+1]; 
end
if nargin < 3 || isempty(lambdapix)
    lambdapix = [805 imageSize(2)/2+1]; 
end
if nargin < 4 || isempty(MmPerPixel)
    MmPerPixel = 0.0104;% * binning; %measured w scale 27/1/25 from getMmPerPix.m
end
if nargin < 5 || isempty(binning)
    binning = 1;
end
if nargin < 6
    topLeftPix = [2048/2/binning - imageSize(1)/2+1 2048/2/binning - imageSize(2)/2+1];
    % [260 188]
end

info.binning = binning;
info.imageSize = imageSize;
info.bregmapix = bregmapix;
info.lambdapix = lambdapix;
info.MmPerPixel = MmPerPixel;
info.topLeftPix = topLeftPix;


