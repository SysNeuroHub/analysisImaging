function image2DMDpng(thisImage, filename)
%create a png while preserving the original pixel size
%with gray colorscale

[height, width] = size(thisImage);

h=figure('visible','off');
h.InnerPosition = [1 1 width height];
imagesc(thisImage);colormap(gray);

axis image off
xlim([1 size(thisImage,2)]);
ylim([1 size(thisImage,1)]);

ax = gca;
ax.Position = [0 0 1 1];

%screen2png(pngName,h);
set(h,'InvertHardcopy','off'); %31/1/22
set(0, 'currentfigure', h); %28/3/18

set(h,'Units','pixels');
scrpos = get(h,'Position');
newpos = scrpos/100;
set(h,'PaperUnits','inches',...
'PaperPosition',newpos)
print('-dpng', filename, '-r100'); %'-r300' will yield higher resolution than the input
%to save as bmp (but the resulting file is heavier than png):
% F = getframe(h); imwrite(F.cdata, [filename '.bmp']);

close(h);
