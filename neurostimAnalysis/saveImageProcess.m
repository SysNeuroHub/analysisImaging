function [imageProc, cic, stimInfo] = saveImageProcess(expInfo, rescaleFac, ...
    rebuildImageData, doRegistration)
% imageProc = saveImageProcess(expInfo, procParam, rebuildImageData,...
%     makeMask, uploadResult)
% save imageProc struct, after temporally aligining oephys & neurostim events
% imageProc contains V that is pixels x frames
% imageProc.V
% imageProc.OETimes
% imageProc.cic
% imageProc.stimInfo
% imageProc.rescaleFac
% imageProc.OEInfo
% imageProc.nrRepeats
% imageProc.conditions
%
% 2022-02-09 DS from quickAnalysisCJ231
% 2022-04-22 temporal filtering is no longer done in this function (to
% reduce computation time)
%
% TODO: replace save names with those from getDataPaths.m

showFig = 0;

addpath('\\storage.erc.monash.edu\shares\R-MNHS-Syncitium\Shared\Daisuke\sandbox');
if exist('C:\Users\dshi0006\git','dir')
    addpath(genpath('C:\Users\dshi0006\git'));
else
    addpath(genpath('C:\git'));    %dsbox,oi-tools,analysis-tools,analysisImaging, npy-matlab,neurostim, marmolab-stimuli
end

if isempty(rebuildImageData)
    rebuildImageData = false;
end

%rescaleFac = procParam.rescaleFac;
%cutoffFreq = procParam.cutoffFreq;
%lpFreq = procParam.lpFreq;

%addDirPrefs;
dirPref = getpref('nsAnalysis','dirPref');
nsOriServer = dirPref.nsOriServer;%\\storage.erc.monash.edu\shares\R-MNHS-Syncitium\Shared\MarmosetData';
oiOriServer = dirPref.oiOriServer;
rootDir = dirPref.rootDir;
oeOriServer = dirPref.oeOriServer; %direct ethernet connection
        
subject = expInfo.subject;
YYYYMMDD = expInfo.date;
nsName = expInfo.nsName;
expID = expInfo.expID;

%% DO NOT CHANGE BELOW
%% oephys digital only
OEInfo.trCh = 1;
OEInfo.camStrobeCh = 7;
OEInfo.ventilatorCh = 4;
OEInfo.expCh = 5;

%% oephys analog
OEInfo.pdCh = 1;

%file name notation
%imaging: ['exp' (consecutive experiment number of an animal)];

%% fix OE name
%fixoephys(fullfile(nsOriServer,expDate)); this is a function to fix OE directory not NS filename

%% rename OEphys filenames (NS data is kept intact)
expDate = sprintf('%s\\%s\\%s',YYYYMMDD(1:4),YYYYMMDD(5:6), YYYYMMDD(7:8));
expDate  = strrep(expDate,'\',filesep);
%nsOriDir = fullfile(nsOriServer, expDate);
oeOriDir = fullfile(oeOriServer, expDate);
fixoephys(oeOriDir,'verbose',true);%skip already renamed folders

%% copy oephys/ns/imaging data to local
DirBase = fullfile(rootDir,expDate);
if ~exist(DirBase,'dir') mkdir(DirBase); end;
cd(DirBase);
fullOEName = retrieveFullOEName(oeOriDir, nsName);


expName = num2str(expID);
OIName = ['exp' num2str(expName)];


%% copy ns data to local
copyServer2Local(nsOriServer, rootDir, fullOEName);

%% copy oephys data to local
copyServer2Local(oeOriServer, rootDir, fullOEName);

%% copy imaging data to local
copyServer2Local(oiOriServer, rootDir, fullOEName, OIName);

%% at this moment, assume all data is uploaded to "rootDir"

DirBase = fullfile(rootDir,expDate);%,expName);
ephysDirBase = strrep(fullfile(DirBase,fullOEName,'Record Node 103\experiment1\recording1'),'\',filesep);
saveDirBase = fullfile(rootDir,subject,'processed');
OEInfo.jsonFile = fullfile(ephysDirBase,'structure.oebin');
stimName = [fullOEName(1:end-20) '.mat'];
stimFile = fullfile(DirBase, stimName);

%% check if these data files exists
saveDir_full = saveDirBase;
imagingDir_full = fullfile(DirBase, OIName);
imDataName = ['imageData_' regexprep(expDate,filesep,'_') '_' expName '_resize' num2str(rescaleFac*100) '.mat'];
imageSaveName = fullfile(saveDirBase, imDataName);

if ~exist(imagingDir_full, 'dir')
    error(['DONT EXIST imaging ' imagingDir_full]);
end
if ~exist(OEInfo.jsonFile, 'file')
    error(['DONT EXIST ephys ' OEInfo.jsonFile]);
end
if ~exist(stimFile, 'file')
    error(['DONT EXIST stim ' stimFile]);
end


%% load camera data
if exist(imageSaveName,'file') && rebuildImageData
    delete(imageSaveName);
end
if exist(imageSaveName,'file')
    disp(['Loading ' imageSaveName]);
    load(imageSaveName,'imageData');
else
    
    imageData = buildImageDataOI(imagingDir_full, rescaleFac, 0, 0, false, ...
        doRegistration);
    
    mkdir(fileparts(imageSaveName));
    save(imageSaveName, 'imageData', '-v7.3');
end

%% load stimulus data
load(stimFile,'c');
stimInfo = getStimInfo(c);
nrRepeats = c.nrTrials;%    c.blocks.nrRepeats;
conditions = get(c.prms.condition, 'atTrialTime', inf);


%% retrieve time stamps
disp('Retrieving timestamps in OpenEphys')
PDTh = []; %CJ246
OETimes = getOETimes(OEInfo, nrRepeats, showFig,PDTh);
close all;

if nrRepeats > numel(OETimes.stimOnTimes)
    disp('Reducing stimInfo.stimLabels');
    stimInfo.stimLabels = stimInfo.stimLabels(1:numel(OETimes.stimOnTimes));
end

%% align timestamps between camera and oephys ... always happens if camera is operated through network
if length(OETimes.camOnTimes) < size(imageData.imstack,3)
    disp('WARNING: IMSTACK LENGTH ADJUSTED');
    %extremely ugly hack to align timestamps between camera and blackrock
    imageData.imstack = imageData.imstack(:,:,1:length(OETimes.camOnTimes));
    save(imageSaveName, 'imageData','-v7.3');
elseif length(OETimes.camOnTimes) > size(imageData.imstack,3)
    disp('WARNING: OETimes.camOnTimes ADJUSTED');
    OETimes.camOnTimes = OETimes.camOnTimes(1:size(imageData.imstack,3));
end
nFrames = length(OETimes.camOnTimes);
Fcam = 1/median(diff(OETimes.camOnTimes)); %camera effective sampling rate


[pspecBefore, axisPspec] = pmtm(imageData.imageMeans(1:nFrames) - mean(imageData.imageMeans(1:nFrames)), ...
    3,nFrames,Fcam);

V = reshape(imageData.imstack, size(imageData.imstack,1)*size(imageData.imstack,2), []);
imageData.imstack = [];

% imageMeansAfter = mean(V,1);
% [pspecAfter, axisPspec] = pmtm(imageMeansAfter-mean(imageMeansAfter), ...
%     3,nFrames,Fcam);


%% mean image across time
figure;
imagesc(imageData.meanImage);
%hold on
%contour(nanMask==1,1,'color','r');
colormap(gray)
axis equal tight
mcolorbar;
saveas(gcf, fullfile(saveDir_full,['meanImage' stimName(1:end-4) '.png']));
close;

%% image means
figure('position',[0 0 1920 1080]);
ax(1)=subplot(211);
plot(OETimes.camOnTimes, imageData.imageMeans(1:nFrames));hold on
% plot(OETimes.camOnTimes, imageMeansAfter);

vbox(OETimes.stimOnTimes, OETimes.stimOffTimes,[],@cool,conditions);
xlabel('time [s]');
ylabel('imagemeans');
mm=mean(imageData.imageMeans);sd=std(imageData.imageMeans);
title(['mean: ' num2str(mm) ', sd: ' num2str(sd) ' (sd/mean: ' num2str(sd/mm*1e2) ' %)']);

ax(3) = subplot(212);
loglog(axisPspec, pspecBefore);hold on;
% loglog(axisPspec, pspecAfter);
if isfield(OETimes,'ventOnTimes')
    respFreq = 1/median(diff(OETimes.ventOnTimes));
    vline(respFreq,ax(3),'-',[.8 .8 .8]);
end
xlabel('Frequency [Hz]');
ylabel('PSD');
%legend('before filtering','after filtering');

saveas(gcf, fullfile(saveDir_full,['imageMeans' stimName(1:end-4) '.png']));
close;

%% save result
imageProc.V = V;
imageProc.OETimes = OETimes;
cic = c;
imageProc.cic = cic;
imageProc.stimInfo = stimInfo;
imageProc.rescaleFac = rescaleFac;
imageProc.OEInfo = OEInfo;
%imageProc.nanMask = nanMask;
%imageProc.procParam = procParam;
imageProc.nrRepeats = nrRepeats;
imageProc.conditions = conditions;

save(imageSaveName,'imageProc','-append');
%save(imageSaveName,'imageProc','cic','stimInfo','-append');


% [~, winSamps, singlePeriEventV, stimLabels_ok, uniqueLabels] ...
%     = eventLockedAvg(V, OETimes.camOnTimes, OETimes.stimOnTimes, ...
%     stimInfo.stimLabels, calcWin);
% %< only use events whose timeWindow is within the recording
% 
% singlePeriEventStack = shiftdim(reshape(singlePeriEventV, numel(stimLabels_ok), ...
%     imageSize_r(1),imageSize_r(2),[]),1);
% %singlePeriEventStack: x - y - t - events
% clear singlePeriEventV
% 
% %% compute dI/I
% baseWin = ones(1,numel(winSamps));
% singlePeriEventStack = getdFFSingleEventStack(singlePeriEventStack, baseWin,...
%     dFFmethod, doMedian);
% 
% %% save stackset 16/11/22
% stackServer = 'E:';
% savgName = ['Stack_' regexprep(expDate, '\','_') '_' expName ...
%     '_resize' num2str(rescaleFac*100) '.mat'];
% 
% tools.SaveMyStacks(singlePeriEventStack, stackServer,'E:\tmp',savgName, 1);


