function image_resized_cropped = resizeCropFullImg(image_full, camImg)
%image_resized_cropped = resizeCropFullImg(image_full, camImg)
%returns image resized and cropped according to camImg
%image_full:2048 x 2048 pixels without binning

resized = imresize(image_full, 1/camImg.binning);

rowIdx = camImg.topLeftPix(1):camImg.topLeftPix(1)+camImg.imageSize(1)-1;
colIdx = camImg.topLeftPix(2):camImg.topLeftPix(2)+camImg.imageSize(2)-1;

image_resized_cropped = resized(rowIdx, colIdx);
