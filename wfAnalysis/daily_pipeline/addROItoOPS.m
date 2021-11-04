function ops = addROItoOPS(ops)
%draw roi as a 2d binary image and append it to ops
%2/11/2021 created from runSVDKT.m

theseFiles = generateFileList(ops, 1);

switch ops.rawDataType
    case 'tif'
        imageForROI = mean(loadTiffStack(theseFiles{1}, 'tiffobj',0),3); %make purple and blue image separately??
    case 'OI'
        imageForROI = mean(loadBlk(theseFiles{1}),3);
end
imagesc(imageForROI);
axis equal tight;
colormap(gray);
caxis(prctile(imageForROI(:),[1 99]));
%colormouse;%change the colormap range interactively using the mouse
title('left click to put anchor points., double click to finish');
roiAhand = images.roi.AssistedFreehand;
draw(roiAhand);
ops.roi = createMask(roiAhand);
close;