function detrendDatFile(datPath, outFile, param2D, imageSize, nFr, ops)

% created from registerDatFile
format = '*uint16'; %28/1/25

batchSize = ops.nRegisterBatchLimit; % images at once
numBatches = ceil(nFr/batchSize);
fid = fopen(datPath);
fidOut = fopen(outFile, 'w');

param1D.k = param2D.k(:);
param1D.beta = param2D.beta(:);
param1D.scaleFac = param2D.scaleFac(:);
try
    for b = 1:numBatches
        imstack = fread(fid,  imageSize(1)*imageSize(2)*batchSize, format);
        imstack = single(imstack);
        imstack = reshape(imstack, imageSize(1)*imageSize(2), []);
        frameNumbers = (1+(b-1)*batchSize):min(b*batchSize, nFr);
        if numel(frameNumbers)>size(imstack,2)
            frameNumbers = frameNumbers(1:size(imstack,2));
        end
        detrendFrames = applyStretchedExp(imstack, param1D, frameNumbers);
        detrendFrames = reshape(detrendFrames, imageSize(1), imageSize(2), []);

        fwrite(fidOut, uint16(detrendFrames), 'uint16');
    end
catch me
    
    fclose(fid);
    fclose(fidOut);
    
    rethrow(me);
end
fclose(fid);
fclose(fidOut);
    
    