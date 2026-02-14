function saveEveryImages(imageStack, saveDir, binary, imgSuffix)
%saveEveryImages(imageStack, saveDir)
%saves every images in imageStack as individual png files
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
    
    % axis ij image off
    % xlim([1 width]);
    % ylim([1 height]);
    % 
    % axpatch = gca;
    % axpatch.Position = [0 0 1 1];
    % 
     % thisName = [num2str(ii)];
    % screen2png(fullfile(saveDir, thisName),f);
    close(f);

end

set(0, 'DefaultFigureVisible', 'on')