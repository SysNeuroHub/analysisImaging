function dat2mov(datPath, outFile, imageSize, nFr, batchSize, crange)
% make .avi uncompressed movie from .dat file
% 26/5/21 created from registerDatFile

if nargin < 6
    crange = [];
end
if nargin < 5
    batchSize = 750;
end
%batchSize = ops.nRegisterBatchLimit; % images at once
numBatches = ceil(nFr/batchSize);
fid = fopen(datPath);
% fidOut = fopen(outFile, 'w');
v = VideoWriter(outFile,'Uncompressed AVI');
open(v);
try
    for b = 1:numBatches
        
        fprintf('%d / %d \n',b, numBatches);
        imstack = fread(fid,  imageSize(1)*imageSize(2)*batchSize, '*int16');
        imstack = single(imstack);
        imstack = reshape(imstack, imageSize(1), imageSize(2), []);
        %         regFrames = register_movie(imstack, ops, ds);
        %         fwrite(fidOut, int16(regFrames), 'int16');
    
        if b==1 && isempty(crange)
            maxval = max(imstack(:));
            minval = 0;
        else
            maxval = crange(2);
            minval = crange(1);
        end
        
        imstack(imstack>maxval) = maxval;
        imstack(imstack<minval) = minval;
        
        imstack = (imstack - minval) / (maxval - minval);
%         imstack = imstack/maxval;
         imstack(imstack>1) = 1;
         imstack(imstack<0) = 0;
        
        for tt = 1:size(imstack,3)
            writeVideo(v, imstack(:,:,tt));
        end
    end
catch me
    
    fclose(fid);
    %fclose(fidOut);
    close(v)
    rethrow(me);
end
fclose(fid);
%fclose(fidOut);
close(v)
    
    