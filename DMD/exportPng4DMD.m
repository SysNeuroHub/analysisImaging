function exportPng4DMD(saveName, fig, binary)
% exportPng4DMD(saveName, fig, binary)
% exports a figure with the original image resolution in fig
% created from showAllenCCFBregmaLambda
% if binary, the resulting image would be 0/255 so png image is readily visible

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
C(1,:) = 0; C(end,:) = 0; %HACK
if binary
    imwrite(intmax("uint8")*uint8(imbinarize(sum(C,3))),[saveName '.png'],'png');
else
    imwrite(im2gray(C),[saveName '.png'],'png');
end
