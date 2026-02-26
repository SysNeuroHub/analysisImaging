function saveEveryImages(imageStack, saveDir, binary, imgSuffix)
%saveEveryImages(imageStack, saveDir)
%saves every images in imageStack as individual png files
% if binary and imageStack is [0 1], it will be expanded to [0 255]
%
% TODO: save as a tiff stack?
% 4/2/26 now use exportPng4DMD

if nargin < 4
    [~,dirname] = fileparts(saveDir);
    imgSuffix = dirname;
end
if nargin < 3
    binary = 0;
end

if binary && max(imageStack(:)) <= 1
    imageStack = uint8(round(double(intmax("uint8"))*imageStack));
end

set(0, 'DefaultFigureVisible', 'off');

width = size(imageStack,2);
height = size(imageStack,1);
nImages = size(imageStack, 3);
for ii = 1:nImages
    disp([num2str(ii) '/' num2str(nImages)]);
    thisName = [num2str(ii) '_' imgSuffix];
    
    f=figure;
    f.InnerPosition = [1 1  width height];

    image(squeeze(imageStack(:,:,ii)));colormap(gray);

    exportPng4DMD(fullfile(saveDir, thisName), f, binary);
    close(f);
end

set(0, 'DefaultFigureVisible', 'on')