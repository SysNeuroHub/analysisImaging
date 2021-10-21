function theseFiles = generateFileList(ops, vidNum)
% theseFiles = generateFileList(ops, vidNum)
% creates a list of video files (tif/mat/bin) under ops.fileBase/ops.inclExpList
%
% ops contains:
%   theseFiles(optional): if exists, use these fullpaths
%   inclExpList: experiment number(s) used to create the list
% 
% called in pipelineHereKT.m


if nargin < 2
    vidNum = 1;
end

theseFiles = {};

if isfield(ops, 'theseFiles') && ~isempty(ops.theseFiles)
    % user has provided a specific list of files to use
    theseFiles = ops.theseFiles;
    
elseif isfield(ops, 'inclExpList') && ~isempty(ops.inclExpList)
    % user has specified specific experiments to include. In this case
    % we're going to look for files in subdirectories corresponding to
    % those experiments
    
    if ops.verbose
        fprintf(ops.statusDestination, 'found subdirectories in fileBase, will look for files there\n');
    end
    for subd = ops.inclExpList
        subFileBase = fullfile(ops.fileBase, num2str(subd));
        if exist(subFileBase)
            subDTheseFiles = directoryFileList(subFileBase, ops.rawDataType);
            theseFiles = [theseFiles subDTheseFiles];
        else
            fprintf(ops.statusDestination, 'warning! Looked for files in %s but that directory did not exist.\n', subFileBase);
        end
    end
    
    
else
    
    % first check if there are subdirectories here. If so we're going to
    % use all the images in them. If not we use the images that are in this
    % folder itself.
    d = dir(fullfile(ops.vids(vidNum).fileBase)); 
    if length(d)>2
        d = d(3:end); % first two are . and ..
        isdir = [d.isdir];
        if sum(isdir)>0
            % has subdirectories
            if ops.verbose 
                fprintf(ops.statusDestination, 'found subdirectories in fileBase, will look for files there\n');
            end
            
            for subd = 1:sum(isdir)
                subFileBase = fullfile(ops.fileBase, d(subd).name);
                subDTheseFiles = directoryFileList(subFileBase, ops.rawDataType);
                theseFiles = [theseFiles subDTheseFiles];
            end
        else
            % no subdirectories, so use files that are here directly
            theseFiles = directoryFileList(ops.vids(vidNum).fileBase, ops.rawDataType);
        end
    end
end


            
function theseFiles = directoryFileList(fileBase, rawDataType)
                
switch rawDataType
    case 'tif'
        theseFilesDir = dir(fullfile(fileBase, '*.tif'));
        [~,ii] = sort([theseFilesDir.datenum]);
        theseFilesDir = theseFilesDir(ii);
    case 'customPCO'
        theseFilesDir = dir(fullfile(fileBase, '*.mat'));
        [~,ii] = sort([theseFilesDir.datenum]);
        theseFilesDir = theseFilesDir(ii);
    case 'StackSet'
        theseFilesDir = dir(fullfile(fileBase, '*.bin'));
    case 'OI'
        theseFilesDir = dir(fullfile(fileBase, '*.BLK'));
        
        % extract expPrefix, expId and blkId from each filename
        % taken from blkImport.m
        finfo = arrayfun(@(x) regexp(x.name,'(?<expPrefix>.*)_E(?<expId>\d+)B(?<blkId>\d+)\.BLK','tokens'), theseFilesDir);

        for expPrefix = unique(cellfun(@(x) x(1), finfo))',
            i = find(cellfun(@(x) strcmp(x(1),expPrefix),finfo));
            
            expIds = cellfun(@(x) str2num(x{2}), finfo(i));
            blkIds = cellfun(@(x) str2num(x{3}), finfo(i));

            orderedIdx = zeros(length(i),1);
            for expId = unique(expIds)',
                j = find(expIds == expId);
                
                [~,k] = sort(blkIds(j)); % try and preserve order of acquisition...
                
                orderedIdx(j) = j(k);
            end
        end
        
    theseFilesDir = theseFilesDir(orderedIdx);
        
end
theseFiles = cellfun(@(x)fullfile(fileBase,x),{theseFilesDir.name},'UniformOutput', false);