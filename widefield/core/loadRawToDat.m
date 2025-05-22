function dataSummary = loadRawToDat(ops)
% dataSummary = loadRawToDat(ops)
% converts a set of tif files in a directory (specified by ops.fileBase) to a
% flat binary (dat) file in datPath. While doing so, extract the data from
% the image stamps about time and frame number, and compute the mean of
% each image and also the mean image. Bin data if requested.
%
% Importantly, this function will try to identify locations of different
% recordings within the files. It does this by assuming any 2-second long
% gap is a break between recordings (the logic is: missed frames appear to
% be rare and just a few at a time so you shouldn't ever miss 2 seconds of
% frames in a row. On the other hand, timeline start/stop always takes at
% least two seconds so you will have at least that much time between
% recordings).
%
% this function assumes every experiment to start with the same LED
%
% ops is a struct that includes:
% - theseFiles - a cell array of filenames
% - datPath - a filename where you would like to create the raw binary file
% - verbose - logical, whether to display status messages
% - rawDataType - can be "tif" or "customPCO", specifies how to load the
% images
% - hasBinaryStamp - specifies whether the raw images have data encoded in
% binary stamps (these are PCO Edge cameras)
% - hasASCIIstamp - specifies whether there are also ASCII stamps which
% will be removed before binning/saving
% - binning - an integer N that specifies to apply N x N binning first
% - frameMod - specifying which images to pull out of the tif files
%    - images included are those for which: mod(frameNumber,
%    frameMod(1))==frameMod(2). So [1 0] takes all images.
% - flipudVid - logical, whether to flip up/down.
extraVerbose = true;

loadMethod = 'tiffobj'; %or 'imread'. NS recommends tiffobj

theseFiles = ops.theseFiles;
datPath = ops.datPath;

frameNumbersFromStamp = [];
timeStampsFromStamp = [];
frameNumbersWithinRec = [];
frameFileIndex = [];
frameRecIndex = [];
imageMeans = [];

lastFilePrefix = []; %8/5/20
lastFileSuffix = []; %8/5/20

recIndex = 0;
lastFrameNum = 0;
total_nfr = 0;
try
    fid = fopen(datPath, 'w');
    frameIndex = 0;
    for fileInd = 1:length(theseFiles)
        
        thisFile = theseFiles{fileInd};
        %[thisFilePrefix, thisFileSuffix] = sscanf(thisFile,'%*s_%d');
        [filepath,fname,fext] = fileparts(thisFile);
        %texts = textscan(fname,'%s %d', 'delimiter', '_'); %for PCO edge
        texts = textscan(fname,'%s %d', 'delimiter', '@'); %for PCO panda
        thisFilePrefix = texts{1}{1};
        thisFileSuffix = texts{2};
        
        if ops.verbose
            fprintf(1, 'loading file: %s (%d of %d)\n', thisFile, fileInd, length(theseFiles));
        end
        
        clear imstack
        switch ops.rawDataType
            case 'tif'
                imstack = loadTiffStack(thisFile, loadMethod, extraVerbose);
            case 'customPCO'
                [~,~,~,imstack] = LoadCustomPCO(thisFile, false, true);
            case 'OI'
                imstack = loadBlk(thisFile);
        end
        
        %if ops.hasBinaryStamp
        if ops.verbose && extraVerbose
            fprintf(1, '  computing timestamps\n');
        end
        
        switch ops.rawDataType
            case {'tif', 'OI'} %rewritten for thorcam
                %records = squeeze(imstack(1,1:14,:));
                %[thisFN, thisTS] = timeFromPCOBinaryMulti(records);
                %thisTS = thisTS*24*3600; % convert to seconds from days
                
                thisFN = 1:size(imstack,3); % + lastFrameNum??
                thisTS =  thisFN + total_nfr; %assuming 1Hz frame rate, i.e. just telling the order
                
                %                     %artificially add 2s gap telling a new recording
                %                      if (~strcmp(thisFilePrefix, lastFilePrefix)) && isempty(thisFileSuffix)
                %                        thisTS = thisTS + 2;
                %                      end
                
                clear theseFrameNumbersWithinRec;
                if fileInd==1
                    firstTS=thisTS(1);
                    thisTS = thisTS-firstTS;
                    %frameDiffs = diff([-10 thisTS]);  % so first frame is a "new rec"
                else
                    thisTS = thisTS-firstTS;
                    %frameDiffs = diff([timeStampsFromStamp(end) thisTS]);
                end
                
                
                
                %newRecInds = find(frameDiffs>=2); % any 2-second gaps between frames indicate a new recording
                if (~strcmp(thisFilePrefix, lastFilePrefix)) && isempty(thisFileSuffix)
                    newRecInds = 1;
                else
                    newRecInds = [];
                end
                theseRecIndex = recIndex*ones(size(thisTS));       %default unless new recs present
                theseFrameNumbersWithinRec = lastFrameNum+(1:length(thisTS)); %default unless new recs present
                nFrThisFile = length(thisTS);
                for n = 1:length(newRecInds)
                    recIndex = recIndex+1;
                    
                    theseFrameNumbersWithinRec(newRecInds(n):nFrThisFile) = 1:(nFrThisFile-newRecInds(n)+1);
                    theseRecIndex(newRecInds(n):nFrThisFile) = recIndex;
                    
                end
                
                inclFrames = mod(theseFrameNumbersWithinRec, ops.frameMod(1))==ops.frameMod(2);
                
                nfr = sum(inclFrames);
                %all consecutive # frames before splitting
                total_nfr = total_nfr + length(theseFrameNumbersWithinRec); %8/5/20
                
                frameNumbersFromStamp(frameIndex+1:frameIndex+nfr) = thisFN(inclFrames); %consecutive frame number across files
                timeStampsFromStamp(frameIndex+1:frameIndex+nfr) = thisTS(inclFrames); %consecutive time stamp across files
                frameFileIndex(frameIndex+1:frameIndex+nfr) = fileInd; %file number
                frameRecIndex(frameIndex+1:frameIndex+nfr) = theseRecIndex(inclFrames); %exp number
                frameNumbersWithinRec(frameIndex+1:frameIndex+nfr) = theseFrameNumbersWithinRec(inclFrames); %frame number within one exp
                lastFrameNum = theseFrameNumbersWithinRec(end);
                
                if ops.frameMod(1)>1
                    if ops.verbose && extraVerbose
                        fprintf(1, '  selecting correct frames\n');
                    end
                    
                end
                
            case 'customPCO'
                fprintf(1, ' customPCO loading not supported in this version. Revert to older code or implement yourself.\n');
                return;
                %                     frameNumbers(frameIndex+1:frameIndex+nfr) = NaN;
                %                     [~,~,thisTS] = LoadCustomPCO(thisFile, false, true);
                %                     thisTS = thisTS*24*3600; % convert to seconds from days
                %
                %                     if fileInd==1
                %                         firstTS=thisTS(1);
                %                     end
                %
                %                     timeStamps(frameIndex+1:frameIndex+nfr) = thisTS-firstTS;
        end
        %         else
        %             %fprintf(1, '  options are set as though you don''t have binary stamps in your images. But you really need them!!! everything is going to fail here.');
        %             %< dont worry Nick, Thorcams can split files of every experiment, that can be used
        %
        %         end
        
        
        
        imstack = removeStamps(imstack, ops.hasASCIIstamp, ops.hasBinaryStamp);
        
        if ops.binning>1
            if ops.verbose && extraVerbose
                fprintf(1, '  binning image\n');
            end
            imstack = binImage(imstack, ops.binning);
        end
        
        if fileInd==1
            sz = size(imstack);
            imageSize = sz(1:2);
            sumImage = zeros(imageSize);
        end
        
        if ops.flipudVid
            imstack = imstack(end:-1:1,:,inclFrames);
        else
            imstack = imstack(:,:,inclFrames);
        end
        
        if ops.verbose && extraVerbose
            fprintf(1, '  computing image means\n');
        end
        imageMeans(frameIndex+1:frameIndex+nfr) = squeeze(mean(mean(imstack,1),2));
        
        
        if ops.verbose && extraVerbose
            fprintf(1, '  computing mean image\n');
        end
        %         sumImage = sumImage+sum(double(imstack),3);
         sumImage = sumImage+sum(imstack,3); %9/1/25
         
        if ops.verbose
            fprintf(1, '  saving to dat\n');
        end
        
        fwrite(fid, imstack, class(imstack));
        
        frameIndex = frameIndex+nfr;
        lastFilePrefix = thisFilePrefix;
        lastFileSuffix = thisFileSuffix;
        
    end
    
catch me
    fclose(fid);
    rethrow(me)
end
fclose(fid);


nFrames = numel(frameNumbersFromStamp);
meanImage = sumImage/nFrames;

dataSummary.frameNumbersFromStamp = frameNumbersFromStamp;
dataSummary.timeStampsFromStamp = timeStampsFromStamp;
dataSummary.frameFileIndex = frameFileIndex;
dataSummary.frameRecIndex = frameRecIndex;
dataSummary.frameNumbersWithinRec = frameNumbersWithinRec;
dataSummary.imageMeans = imageMeans;
dataSummary.meanImage = meanImage;
dataSummary.imageSize = imageSize;
dataSummary.dataType = class(imstack);
dataSummary.nFrames = nFrames;

if ops.verbose
    fprintf(1, '  done, found and loaded %d images from %d recordings\n', nFrames, length(unique(dataSummary.frameRecIndex)));
end
