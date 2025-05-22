function allDs = alignToTarget(datPath, targetFrame, imageSize, nFr, ops)
%allDs = alignToTarget(datPath, targetFrame, imageSize, nFr, ops)
%reads .dat file, figures out the shifts required to align to targetFrame
%
% Inputs:
%   datpath: fullpath to the movie file converted to .dat format
%   targetFrame: image to which image shifts are computed
%   imageSize: image size stored in .dat [pixY pixX]
%   nFr: total number of frames stored in .dat
%   ops: must include 
%       nRegisterBatchLimit: number of frames to process at once
%
% Outputs:
%   allDs: shift size of all frames stored in .dat
%
% called in runSVDKT.m

% format = '*int16';
format = '*uint16'; %28/1/25

batchSize = ops.nRegisterBatchLimit; % images at once
numBatches = ceil(nFr/batchSize);
fid = fopen(datPath);
ind = 1;
allDs = zeros(nFr, 2);
try
    for b = 1:numBatches
        imstack = fread(fid,  imageSize(1)*imageSize(2)*batchSize, format);
        imstack = single(imstack);
        imstack = reshape(imstack, imageSize(1), imageSize(2), []);
        [ds, ~]  = registration_offsets(imstack, ops, targetFrame, 0);

        allDs(ind:ind+size(imstack,3)-1,:) = ds;
        ind = ind+size(imstack,3);
    end
catch me
    
    fclose(fid);
    rethrow(me);
end

fclose(fid);    
    