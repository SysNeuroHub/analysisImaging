

%input 2pStackDir
twopStackDir = '\\vault-v2.erc.monash.edu.au\MNHS-dshi0006\Subjects\TIGRE2GCaMP6s_302\2020-11-25\anterior zstack  607 800';
% wfName = '\\vault-v2.erc.monash.edu.au\MNHS-dshi0006\Subjects\TIGRE2GCaMP6s_302\2020-11-24\05x\blue1360ms.tif';
wfName = '\\vault-v2.erc.monash.edu.au\MNHS-dshi0006\Subjects\TIGRE2GCaMP6s_302\2020-11-24\2x\blue anterior.tif';
rot270 = 1; %rotate wf image 270deg

%% load 2p zstack

% taken from generateFileList/directoryFileList(ops.vids.fileBase, ops.rawDataType);
theseFilesDir = dir(fullfile(twopStackDir, '*.tif'));
[~,ii] = sort([theseFilesDir.datenum]);
theseFilesDir = theseFilesDir(ii);
these2pFiles = cellfun(@(x)fullfile(twopStackDir,x),{theseFilesDir.name},'UniformOutput', false);

% taken from loadRawToDat
stack2p = [];
for fileInd = 1:length(these2pFiles)
    
    thisFile = these2pFiles{fileInd};
    %[thisFilePrefix, thisFileSuffix] = sscanf(thisFile,'%*s_%d');
    [filepath,fname,fext] = fileparts(thisFile);
    texts = textscan(fname,'%s %d', 'delimiter', '_');
    thisFilePrefix = texts{1}{1};
    thisFileSuffix = texts{2};
    stack2p(:,:,fileInd) = loadTiffStack(thisFile, 'imread', 1);
end
twop = mean(stack2p,3);        



%% load wf image
wf = loadTiffStack(wfName, 'imread',1);
if rot270
    wf = rot90(wf, -1);
end

imagesc(wf)


   