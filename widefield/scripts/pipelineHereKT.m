function pipelineHereKT()
%this script does
%image registration if ops.doRegistration


% try % putting the entire script in a try-catch block to we can return even if it fails


% addpath(genpath('\\zserver.cortexlab.net\Code\Rigging\main'));
% addpath(genpath('\\zserver.cortexlab.net\Code\Rigging\cb-tools')); % for some dat package helper functions
% addpath(genpath('C:\Users\Experiment\Documents\GitHub\npy-matlab'))
% addpath(genpath('C:\Users\Experiment\Documents\GitHub\widefield'))
% addpath(genpath('C:\Users\Experiment\Documents\MATLAB\visbox'));
% addpath(genpath('C:\npy-matlab'));
% addpath('C:\svdinput');

close all;

showFig = 1;
saveDat = 1; %if 1, save dat files into the Vault server, else delete them

%serverDir = '\\ad.monash.edu\home\User006\dshi0006\Documents\tempDataServer'; %7/5/20
%'\\lugaro.cortexlab.net\bigdrive\staging\';

load ops.mat; % this must be present in the current directory
diaryFilename = sprintf('svdLog_%s_%s.txt', ops.mouseName, ops.thisDate);
diary(diaryFilename);
    
ops.localSavePath = pathForThisOS(ops.localSavePath);
for v = 1:length(ops.vids)
    ops.vids(v).fileBase = pathForThisOS(ops.vids(v).fileBase);
end

if ~exist(ops.localSavePath, 'dir')
    mkdir(ops.localSavePath);
end
save(fullfile(ops.localSavePath, 'ops.mat'), 'ops');    

%% load all movies into flat binary files, for each color

%load(fullfile(ops.localSavePath, 'results.mat'), 'results');%TEMP START
% for v = 1:length(ops.vids)
%     
%     clear loadDatOps;
%     
%     ops.theseFiles = [];
%     theseFiles = generateFileList(ops, v);
%     
%     ops.vids(v).theseFiles = theseFiles;
%     loadDatOps.theseFiles = theseFiles;
%         
%     ops.vids(v).thisDatPath = fullfile(ops.localSavePath, ['vid' num2str(v) 'raw.dat']);
% end
%% TEMP END

for v = 1:length(ops.vids)
    
    clear loadDatOps;
    
    ops.theseFiles = [];
    theseFiles = generateFileList(ops, v);
    
    ops.vids(v).theseFiles = theseFiles;
    loadDatOps.theseFiles = theseFiles;
        
    ops.vids(v).thisDatPath = fullfile(ops.localSavePath, ['vid' num2str(v) 'raw.dat']);
    loadDatOps.datPath = ops.vids(v).thisDatPath;    
    loadDatOps.verbose = ops.verbose;
    loadDatOps.rawDataType = ops.rawDataType;
    
    loadDatOps.frameMod = ops.vids(v).frameMod;
    loadDatOps.hasASCIIstamp = ops.hasASCIIstamp;
    loadDatOps.hasBinaryStamp = ops.hasBinaryStamp;
    loadDatOps.binning = ops.binning;
    loadDatOps.flipudVid = ops.vids(v).flipudVid;
    
    dataSummary = loadRawToDat(loadDatOps);

   
    if showFig
        %check how frameNumbersFromStamp & frameNumbersWithinRec are used from
        %here
        %sanity check
        subplot(211);
        plot(dataSummary.frameRecIndex); hold on
        plot(dataSummary.frameFileIndex);
        legend('frameRecIndx','frameFileIndex');
        grid on;
        subplot(212);
        plot(dataSummary.imageMeans);
    end
    

    fn = fieldnames(dataSummary);
    results(v).name = ops.vids(v).name;
    for f = 1:length(fn)
        results(v).(fn{f}) = dataSummary.(fn{f});
    end
    
  %  save(fullfile(ops.localSavePath, 'results.mat'), 'results');
  %  %TEMPORALLY DISABLED
end


%% do image registration? 
% Register a master image and apply the registration to the other movies
if ops.doRegistration

    regOps.NimgFirstRegistration = ops.NimgFirstRegistration;
    regOps.NiterPrealign = ops.NiterPrealign;
    regOps.SubPixel = ops.SubPixel;
    regOps.RegPrecision = ops.RegPrecision;
    regOps.phaseCorrelation = ops.phaseCorrelation;
    regOps.nRegisterBatchLimit = ops.nRegisterBatchLimit;
    regOps.useGPU = ops.useGPU;
    
    v = ops.masterVid;
    datPath = ops.vids(v).thisDatPath;
    
    % determine target frame for the master video
    if ops.verbose
        fprintf(1, 'determining target for image registration\n');
    end
    imageSize = results(v).imageSize;
    nFr = results(v).nFrames;
    targetFrame = determineTargetFrame(datPath, imageSize, nFr, regOps);
    
    
    % figure out the shifts required to align to it
    if ops.verbose
        fprintf(1, 'determining registration shifts\n');
    end
    ds = alignToTarget(datPath, targetFrame, imageSize, nFr, regOps);
    
    % now shift every video to match
    if ops.verbose
        fprintf(1, 'applying registration\n');
    end
    for v = 1:length(ops.vids)
        results(v).registrationDs = ds;
        if ops.verbose
            fprintf(1, '  to vid %d\n', v);
        end
        datPath = ops.vids(v).thisDatPath;
        regPath = fullfile(ops.localSavePath, ['vid' num2str(v) 'reg.dat']);
        ops.vids(v).thisRegPath = regPath;
        registerDatFile(datPath, regPath, ds, results(v).imageSize, results(v).nFrames, regOps);
    end
    
    save(fullfile(ops.localSavePath, 'results.mat'), 'results');%16/10/20

end

%% do hemodynamic correction?
% - don't do this here - it likely works just as well on SVD representation
% (though that has not been explicitly tested). 



%% perform SVD

svdOps.NavgFramesSVD = ops.NavgFramesSVD;
svdOps.verbose = ops.verbose;
svdOps.nSVD = ops.nSVD;
svdOps.useGPU = ops.useGPU;

% If an ROI for the brain was selected to exclude outside pixels
% (AP 160929)
% < better to do this before image registration?
if isfield(ops,'roi')
    svdOps.roi = ops.roi;
end

for v = 1:length(ops.vids)
    fprintf(1, ['svd on ' ops.vids(v).name '\n']);
    
    svdOps.Ly = results(v).imageSize(1); svdOps.Lx = results(v).imageSize(2); % not actually used in SVD function, just locally here

    if ops.doRegistration 
        %minDs = min(results(v).registrationDs, [], 1);
        %maxDs = max(results(v).registrationDs, [], 1);

        %svdOps.yrange = ceil(maxDs(1)):floor(svdOps.Ly+minDs(1)); %correct?
        %svdOps.xrange = ceil(maxDs(2)):floor(svdOps.Lx+minDs(2)); %correct?   
        
        svdOps.yrange = 1:svdOps.Ly; %20/10/20
        svdOps.xrange = 1:svdOps.Lx; %20/10/20
        svdOps.RegFile = ops.vids(v).thisRegPath;
    else
        svdOps.yrange = 1:svdOps.Ly; % subselection/ROI of image to use
        svdOps.xrange = 1:svdOps.Lx;
        svdOps.RegFile = ops.vids(v).thisDatPath;
    end
    svdOps.Nframes = results(v).nFrames; % number of frames in whole movie

    svdOps.mimg = results(v).meanImage;
        
    tic
    [U, Sv, V, totalVar] = get_svdcomps(svdOps);
    
    %     Error using reshape
    %     Product of known dimensions, 250000, not divisible into total number of elements,
    %     245301504.
    %
    %     Error in get_svdcomps (line 68)
    %     data = reshape(data, Ly, Lx, []);

    toc   
    
    % what to do about this? Need to save all "vids" - where?
    fprintf(1, 'attempting to save to server\n');
    
    saveOps.rigName = ops.vids(v).rigName;
    saveOps.verbose = ops.verbose;
    saveOps.mouseName = ops.mouseName;
    saveOps.thisDate = ops.thisDate;
    saveOps.vidName = ops.vids(v).name;
    saveOps.saveAsNPY = ops.saveAsNPY;
    saveOps.frameMod = ops.vids(v).frameMod;
    saveOps.expRefs = ops.expRefs;
    saveOps.inclExpList = ops.inclExpList; %13/5/20
    
    results(v).Sv = Sv;
    results(v).totalVar = totalVar;
    
    %% match timestamps with Timeline, split V into each exp, 
    % and upload to the market server specified as
    % dat.expPath(ops.mouseName, thisDateStr, thisSeriesNum, 1, 'main', 'master')

     V = V(:,1:43518); %hack for hercules exp3

    saveSVD(saveOps, U, V, results(v)); 
    
    results(v).U = U;
    results(v).V = V;
    
    clear U V %9/11/20
    
end

inspectSVDresult(results); %15/5/20
filePath = dat.expPath(saveOps.expRefs{1}, 'main', 'master');
movefile('*.png', filePath);

%% save results.mat locally commented out 3/2/25
% fprintf(1, 'saving all results locally at %s\n', fullfile(ops.localSavePath, 'results.mat'));
% save(fullfile(ops.localSavePath, 'results.mat'), 'results', '-v7.3');
% fprintf(1, 'done saving results.mat\n');

rng('shuffle','twister');

% if isfield(ops, 'emailAddress') && ~isempty(ops.emailAddress)        
% 
%     messages = {'I am the SVD master.', 'Decomposing all day, decomposing all night.', 'Yes! Yes! Woooooooooo!', 'Wha wha whaaat?? It happened! It really happened!!'};    
%     % Send the email
%     mailFromLugaro(ops.emailAddress, [ops.mouseName '_' ops.thisDate ' finished.'], ...
%         messages{randi(numel(messages),1)}, diaryFilename);
% 
% end

% save(fullfile(ops.localSavePath, 'done.mat'), []);
% Instead, copy the folder of raw files into the /mnt/data/toArchive folder

% Diary has to be turned off before moving, otherwise it'll move the diary
% and error out when it tries to write to it
diary off;

% clean up dat files or move them to the Vault serer
for v = 1:length(ops.vids)
    if isfield(ops.vids(v), 'thisDatPath') && exist(ops.vids(v).thisDatPath) ...
            && ~strcmp(ops.vids(v).thisDatPath, ops.vids(v).fileBase)
        if saveDat
            continue;
        else
            delete(ops.vids(v).thisDatPath);
        end
    end
    if isfield(ops.vids(v), 'thisRegPath') && exist(ops.vids(v).thisRegPath) ...
            && ~strcmp(ops.vids(v).thisDatPath, ops.vids(v).fileBase)
        if saveDat
             continue;
        else
            delete(ops.vids(v).thisRegPath);
        end
    end
end

% % %% Copy rawdata and processed data into vault server
% % %destFolder = fullfile(serverDir, [ops.mouseName '_' ops.thisDate]);
% % thisDateStr = ops.thisDate(1:10); %8/5/20
% % thisSeriesNum = str2num(ops.thisDate(12:end)); %8/5/20
% % destFolder = fileparts(dat.expPath(ops.mouseName, thisDateStr, thisSeriesNum,1,'vault','master'));
% % 
% % % if ~strcmp(destFolder, ops.vids(1).fileBase) %commented out 20/1/20
% % fprintf(1, 'copying files to vault server\n');
% % mkdir(destFolder);
% % if  isfield(ops, 'inclExpList') && ~isempty(ops.inclExpList) && ~isequal(ops.vids(1).fileBase, destFolder)
% %     for ee = 1:length(ops.inclExpList)
% %         copyfile(fullfile(ops.vids(1).fileBase, num2str(ops.inclExpList(ee))),...
% %             fullfile(destFolder, num2str(ops.inclExpList(ee)))); %will change to movefile
% %         
% %     end
% %     %movefile can produce error. better use copy here then rmdir later
% %     copyfile(fullfile(ops.localSavePath,'*'), ...
% %         fullfile(destFolder, num2str(ops.inclExpList(ee))));
% %     
% % else
% %     if  ~isequal(ops.vids(1).fileBase, destFolder)
% %         copyfile(fullfile(ops.vids(1).fileBase, '*'), destFolder); %will change to movefile
% %     end
% %     copyfile(fullfile(ops.localSavePath,'*'), destFolder); 
% % end
% % 
% % 
% % % end
% % 
% % if ~isequal(ops.vids(1).fileBase, destFolder)
% %     try
% %         for ee = 1:length(ops.inclExpList)
% %             delete(fullfile(ops.vids(1).fileBase, num2str(ops.inclExpList(ee)),'*')); %25/6/20
% %             disp(['Deleted ' fullfile(ops.vids(1).fileBase, num2str(ops.inclExpList(ee)))]);
% %         end
% %     catch err
% %         disp(['Could not delete '  fullfile(ops.vids(1).fileBase, num2str(ops.inclExpList(ee)))]);
% %     end
% % end
% % fprintf(1, 'done moving files\n');
% % 
% % 
% % %% delete the working folder
% % if length(dir(ops.localSavePath))==2 % if the folder is empty
% %     %this condition ensures that the all the data is moved to the server
% % 
% %     cd(fileparts(ops.localSavePath));
% %     rmdir(ops.localSavePath);
% %     fprintf('moved the process files to the Vault server');
% % end
% % 
% % fprintf(1, 'done backup and cleanup\n');

% catch me
%     disp(me.message);
% %     if isfield(ops, 'emailAddress') && ~isempty(ops.emailAddress)
% %         if exist(diaryFilename)
% %             mailFromLugaro(ops.emailAddress, [ops.mouseName '_' ops.thisDate ' got an error :('], ...
% %                 me.message, diaryFilename);
% %         else
% %             mailFromLugaro(ops.emailAddress, [ops.mouseName '_' ops.thisDate ' got an error :('], ...
% %                 me.message, []);
% %         end
% %     end
%     diary off;
% end
% 
% end