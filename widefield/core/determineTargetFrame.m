function targetFrame = determineTargetFrame(datFile, imageSize, nFr, regOps)
% targetFrame = determineTargetFrame(datFile, imageSize, nFr, regOps)
% returns a average image (after registration) across randomly sampled images in datFile
%
% Inputs:
%   datFile: fullpath to the movie file
%   imageSize: number of pixels housed [Y X] in the datFile
%   nFr: number of frames housed in the datFile
%   regOps: ops passed into align_iterative
%
% Output:
%   targetFrame: average image across randomly sampled images in datfile

imgInds = 1:ceil(nFr/regOps.NimgFirstRegistration):nFr; % want these frames, evenly distributed from the recording

m = memmapfile(datFile, 'Format', {'int16' [imageSize(1) imageSize(2) nFr] 'd'});

[AlignNanThresh, ErrorInitialAlign, dsprealign, targetFrame] = ...
    align_iterative(single(m.Data.d(:,:,imgInds)), regOps);

