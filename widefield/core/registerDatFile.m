function registerDatFile(datPath, outFile, ds, imageSize, nFr, ops)
% registerDatFile(datPath, outFile, ds, imageSize, nFr, ops)
% reads .dat file(datpath), align the data according to ds, and save the aligned
% data into another .dat file as "outFile"

% format = '*int16';
format = '*uint16'; %28/1/25

batchSize = ops.nRegisterBatchLimit; % images at once
numBatches = ceil(nFr/batchSize);
fid = fopen(datPath);
fidOut = fopen(outFile, 'w');

try
    for b = 1:numBatches
        imstack = fread(fid,  imageSize(1)*imageSize(2)*batchSize, format);
        imstack = single(imstack);
        imstack = reshape(imstack, imageSize(1), imageSize(2), []);
        regFrames = register_movie(imstack, ops, ds);

%        fwrite(fidOut, int16(regFrames), 'int16');
        fwrite(fidOut, uint16(regFrames), 'uint16');
    end
catch me
    
    fclose(fid);
    fclose(fidOut);
    
    rethrow(me);
end
fclose(fid);
fclose(fidOut);
    
    