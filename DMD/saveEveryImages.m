function saveEveryImages(imageStack, saveDir)
%saveEveryImages(imageStack, saveDir)
%saves every images in imageStack as individual png files
%
% TODO: save as a tiff stack?

width = size(imageStack,2);
height = size(imageStack,1);
nImages = size(imageStack, 3);
for ii = 1:nImages
    f=figure;
        f.InnerPosition = [1 1  width height];

        image(squeeze(imageStack(:,:,ii)));colormap(gray);
    axis ij image off
        xlim([1 width]);
        ylim([1 height]);
        
        axpatch = gca;
        axpatch.Position = [0 0 1 1];

        thisName = [num2str(ii)];
        screen2png(fullfile(saveDir, thisName),f);
        close(f);
end