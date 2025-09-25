<<<<<<< HEAD
function im = plotTVslice(thisSlice)
% function plotTVslice(thisSlice)
% 
% useful way to plot a slice you get from sliceByVector, if the slice was
% made from the allen atlas ccf "template volume"

im = imagesc(thisSlice);
% set(gca, 'YDir','normal')
axis image
axis off
colormap gray
=======
function im = plotTVslice(thisSlice)
% function plotTVslice(thisSlice)
% 
% useful way to plot a slice you get from sliceByVector, if the slice was
% made from the allen atlas ccf "template volume"

im = imagesc(thisSlice);
% set(gca, 'YDir','normal')
axis image
axis off
colormap gray
>>>>>>> origin/main
caxis([0 400]);