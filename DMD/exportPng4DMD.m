function exportPng4DMD(saveName, fig, binary)
% export a figure with the original image resolution in fig
% created from showAllenCCFBregmaLambda
if nargin < 3
    binary = 0;
end
if nargin < 2
    fig = gcf;
end

f=fig;
hImage = findobj(fig, 'Type', 'image');
imageData = get(hImage, 'CData');
[height, width] = size(imageData);

f.InnerPosition = [1 1 width height];
axis ij image off
colormap(gray);
xlim([1 width]);
ylim([1 height]);

ax = gca;
ax.Position = [0 0 1 1];

screen2png(saveName,f); %saved as color


C = imread([saveName '.png']);
if binary
    imwrite(im2bw(C),[saveName '.png'],'png');
else
    imwrite(im2gray(C),[saveName '.png'],'png');
end
