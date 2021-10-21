

function saveSVDtemp(ops, U, V, dataSummary)
% saveSVD(ops, U, V, dataSummary)
% matches timestamps with Timeline, split V into each exp, and upload to
% the market server to: dat.expPath(ops.mouseName, thisDateStr, thisSeriesNum, 1, 'main', 'master')
%
% dataSummary must have:
% - meanImage
% - frameRecIndex, 1xN integers where N is the number of frames in V and
% the entries identify which recording (sequentially) that frame was part
% of


[numExps, nFrPerExp, allT, existExps, alignmentWorked] = determineTimelineAlignmentstemp(ops, size(V,2));

if ops.verbose
    fprintf(1, 'saving SVD results to server... \n');
end

%expNums = dat.parseExpRef(ops.expRefs);%14/5/20

% upload results to server
thisDateStr = ops.thisDate(1:10); %8/5/20
thisSeriesNum = str2num(ops.thisDate(12:end)); %8/5/20
% filePath = dat.expPath(ops.mouseName, thisDateStr, thisSeriesNum, 1, 'main', 'master');
%Upath = fileparts(filePath); % root for the date - we'll put U (etc) and data summary here
filePath = dat.expPath(ops.mouseName, thisDateStr, thisSeriesNum, 4, 'main', 'master'); %6/11/20
Upath = filePath;%6/11/20
% if ~exist(Upath)
%     mkdir(Upath);
% end

if ops.verbose
    fprintf(1, '  saving U... \n');
end

saveU(U, dataSummary.meanImage, Upath, ops);
save(fullfile(Upath, ['dataSummary_' ops.vidName]), 'dataSummary', 'ops');
    

if alignmentWorked

    allV = V;    
    fileInds = cumsum([0 nFrPerExp]);

    for n = 1:numExps
        if ops.verbose
            fprintf(1, '  saving V for exp %s with timestamps... \n', existExps{n});
        end
        filePath = dat.expPath(existExps{n}, 'main', 'master');
        mkdir(filePath);
        
        %V = allV(:,fileInds(n)+1:fileInds(n+1));
        t = allT{n}(1:length(V)); %taken from Timeline via strobeTimes.m

        saveV(V, t, filePath, ops);
        
    end
    % elseif isfield(ops, 'expRefs') && ~isempty(ops.expRefs) && length(ops.expRefs)==length(unique(dataSummary.frameRecIndex))
    %     % here, alignment didn't work but the number of recordings we found
    %     % matches the number of expRefs we were given. So we will save each section of V,
    %     % without timestamps, to its correct directory
    %     allV = V;
    %     recInds = dataSummary.frameRecIndex;
    %     numExps = length(ops.expRefs);
    %
    %     for n = 1:numExps
    %         if ops.verbose
    %             fprintf(1, '  saving V for exp %s without timestamps... \n', ops.expRefs{n});
    %         end
    %         filePath = dat.expPath(ops.expRefs{n}, 'main', 'master');
    %         mkdir(filePath);
    %
    %         V = allV(:,recInds==n);
    %
    %         saveV(V, [], filePath, ops);
    %
    %     end
    
else        
    % alignment didn't work at all, just save it like U, in the root directory
        
    if ops.verbose
        fprintf(1, '  saving V in root... \n');
    end
    
    vPath = Upath;
    
    saveV(V, [], vPath, ops);
    
end
    
% Register results files with database here??

if ops.verbose
    fprintf(1,'done saving\n');
end


function saveU(svdSpatialComponents, meanImage, Upath, ops)

fn = fullfile(Upath, ['svdSpatialComponents_' ops.vidName]);
fnMeanImage = fullfile(Upath, ['meanImage_' ops.vidName]);

if isfield(ops, 'saveAsNPY') && ops.saveAsNPY
    writeUVtoNPY(svdSpatialComponents, [], fn, []);
    writeNPY(meanImage, [fnMeanImage '.npy']);
else
    save(fn, '-v7.3', 'svdSpatialComponents');
    save(fnMeanImage, 'meanImage');
end

function saveV(svdTemporalComponents, t, Vpath, ops)

fn = fullfile(Vpath, ['svdTemporalComponents_' ops.vidName]);
fnT = fullfile(Vpath, ['svdTemporalComponents_' ops.vidName '.timestamps']);

if isfield(ops, 'saveAsNPY') && ops.saveAsNPY
    writeUVtoNPY([], svdTemporalComponents, [], fn);
    if ~isempty(t)
        writeNPY(t, [fnT '.npy']);
    end
else
    if ~isempty(t)
        save(fn, '-v7.3', 'svdTemporalComponents', 't');
    else
        save(fn, '-v7.3', 'svdTemporalComponents');
    end
end
