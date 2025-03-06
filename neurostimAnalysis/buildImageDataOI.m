function imageData = buildImageDataOI(imagingDir_full, rescaleFac, makeROI, ...
    makeAvi, loadFirst, doRegistration)
% imageData = buildImageDataOI(imagingDir_full, makeROI, makeAvi, loadFirst)

%fixed parameters
batchSize = 1e3; %imaging loading

if nargin <2
    rescaleFac = 1;
end
if nargin<3
    makeROI = true;
end
if nargin<4
    makeAvi = false;
end
if nargin<5
    loadFirst = 0;
end
if nargin<6
    doRegistration = false;
end
load('intrinsicImagingOps.mat');
    %where resulting .dat is saved
    ops.localSavePath = fileparts(imagingDir_full); %'E:\tmp';
    
    
    
    if ~exist(ops.localSavePath)
        mkdir(ops.localSavePath);
    end
    
    ops.vids(1).fileBase = imagingDir_full;
    ops.theseFiles = generateFileList(ops);
    disp(ops.theseFiles');
    
    
    [~,expName] = fileparts(imagingDir_full);
    loadDatOps.datPath = fullfile(ops.localSavePath,[expName '_vid' num2str(1) 'raw.dat']);
    if exist(loadDatOps.datPath,'file')
        delete(loadDatOps.datPath);
    end
    if loadFirst
        ops.theseFiles = ops.theseFiles(1);
    end
    loadDatOps.theseFiles = ops.theseFiles;
    loadDatOps.verbose = ops.verbose;
    loadDatOps.rawDataType = ops.rawDataType;
    
    loadDatOps.frameMod = ops.vids(1).frameMod;
    loadDatOps.hasASCIIstamp = ops.hasASCIIstamp;
    loadDatOps.hasBinaryStamp = ops.hasBinaryStamp;
    loadDatOps.binning = ops.binning;
    loadDatOps.flipudVid = ops.vids(1).flipudVid;
    
    if makeROI
        loadDatOps = addROItoOPS(loadDatOps);
    end
    
    %loadRawToDat_ROI in analysisImaging/widefield/core
    imageData = loadRawToDat_ROI(loadDatOps);%, roiy, roix); %only save roi
    
    disp([num2str(imageData.nFrames) ' frames saved to ' loadDatOps.datPath]);
    
    if makeAvi
        %% create video. dat2mov in sandbox
        dat2mov(datPath, [datPath(1:end-3) 'avi'], imageData(1).imageSize, imageData(1).nFrames);
    end
    
    imageSize = imageData(1).imageSize;
    nFr = imageData(1).nFrames;
    numBatches = ceil(nFr/batchSize);
    
    disp(['loading ' loadDatOps.datPath]);
    fid = fopen(loadDatOps.datPath);
    imageData.imstack = [];
    for b = 1:numBatches
        
        fprintf('%d / %d \n',b, numBatches);
        imstack_c = fread(fid,  imageSize(1)*imageSize(2)*batchSize, '*single');%'*int16');
        %imstack_c = single(imstack_c);
        imstack_c = reshape(imstack_c, imageSize(1), imageSize(2), []);
        
        %% image registration
        %\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\master\stacks
        if doRegistration
            if b==1
                [dx, dy, target, imstack_c] = ...
                    tools.RapidReg(imstack_c, 'auto');
            else
                [dx_c, dy_c, target, imstack_c] = ...
                    tools.RapidReg(imstack_c, target);
                dx = [dx_c dx];
                dy = [dy_c dy];
            end
        end
        
        imstack_c = imresize(imstack_c, rescaleFac); %image resize 10/6/21
        
        imageData.imstack = cat(3,imageData.imstack, imstack_c);
    end
    clear imstack_c
    imageSize_r = size(imageData.imstack,[1 2]);
    imageData.imageSize = imageSize_r;
    if doRegistration
        imageData.dx = dx;
        imageData.dy = dy;
    end
    if makeROI
        imageData.mask = imresize(loadDatOps.roi, imageData.imageSize, 'nearest');
        %imageData.imstack = imageData.imstack.*imageData.mask;
    end
    imageData.meanImage = imresize(imageData.meanImage, imageData.imageSize);
    fclose(fid);
    
%     mkdir(fileparts(imageSaveName));
%     save(imageSaveName, 'imageData');
