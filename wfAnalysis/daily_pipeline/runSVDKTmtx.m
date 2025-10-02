%%
% This script requires:
% Timeline data uploaded to master server
%
% This script does:
% 
% apply SVD
% visualize SVD results (inspectSVDresult.m)
% check if #frames are identical that recorded in TL (in saveSVD.m)
% split V into each experiment (in saveSVD.m)
% save U and V to the market server (in saveSVD.m)

% things to do before running this script:
% 1, upload raw data to Market
% 2, download the raw data under E:\Subjects

clear all
addpath 'C:\Documents\git\analysisImaging'
setPath_analysisImaging;

cd('C:\Documents\git\analysisImaging\wfAnalysis\daily_pipeline');
load('amberRedOps.mat');
% For blue/purple alternate
%load('bluePurpleOps.mat'); 
% For purple only
% load('C:\Users\Experiment\Documents\MATLAB\purpleOps.mat')
mouseName = 'test';%'susanoo';
thisDate = '2025-08-14';%'2024-11-22'; %[datestr(now,'yyyy-mm-dd')];  
thisSeries = 1;
expNums = [6];
makeROI = -1;
doRegistration = 0;
hwbinning = 1; %automatically retrieve this from thorcam header??
magnification = .5; 

if ispc
    cd('C:\Documents\git\analysisImaging\wfAnalysis\daily_pipeline');
elseif isunix
    cd('/home/daisuke/Documents/git/analysisImaging/wfAnalysis/daily_pipeline');
end

rootDrive = 'D:\svdinput';
rawDataDir = 'D:\Subjects'; %local temporary storage


%where SVD and summary data is saved
%dat.expFilePath(ops.mouseName, thisDate, thisSeries, expNums, 'widefield','master')

%where the raw and processed data is saved
%dat.expFilePath(ops.mouseName, thisDate, thisSeries, expNums, 'widefield_raw','master')


%% modify ops for this experiment
%hack for 9/6/20 phpL
% ops.vids(1).frameMod = [2 0];
% ops.vids(2).frameMod = [2 1];
%ops.vids = ops.vids(1);%3/7/20

thisDateSeries = [thisDate, '_' num2str(thisSeries)];
if makeROI == 0
    if exist(fullfile(rawDataDir, mouseName,'thisROI.mat'),'file')
        load(fullfile(rawDataDir, mouseName,'thisROI.mat'),'thisROI');
        ops.roi = thisROI;
    else
        ops.roi = [];
    end
end

ops.mouseName = mouseName;
ops.thisDate = thisDateSeries;
%if using inclExpList, specify ops.fileBase
ops.inclExpList = expNums;
ops.fileBase = fullfile(rawDataDir, mouseName, thisDateSeries);
ops.statusDestination = fopen('test.txt','w');
ops.userName = 'Daisuke';
ops = rmfield(ops,'emailAddress');
ops.rigName = 'alloptrig';%'wfrig'; %used in determineTimelineAlignments
ops.doRegistration = doRegistration;
ops.useGPU = 1; %used only for registration?
ops.objectiveType = num2str(magnification); %0.5 / 0.8
ops.pixelSizeUM = 5/str2num(ops.objectiveType)*hwbinning; %added 15/5/20
ops.hasASCIIstamp = 0; %18/5/20


for e = 1:length(expNums)
    expRefs{e} = dat.constructExpRef(mouseName, thisDate, thisSeries, expNums(e));    
end
ops.expRefs = expRefs;

% retrieve camera sampling rate and expos   ure duration from Timeline
timelinePath = dat.expFilePath(expRefs{1}, 'timeline', 'master');
[~,timelineName] = fileparts(timelinePath);
timelinePath = fullfile(ops.fileBase,num2str(expNums(1)),timelineName);
load(timelinePath);
[strobeOnTimes, ~, strobeDurs] = getStrobeTimes(Timeline, ops.rigName);
exposureDur = median(strobeDurs); %[s]
FsPerColor = 1/median(diff(strobeOnTimes))/length(ops.vids);
ops.Fs = FsPerColor;%16/10/20
clear Timeline


for i = 1:length(ops.vids)
    %where video data is stored
    ops.vids(i).fileBase = fullfile(rawDataDir, mouseName, thisDateSeries);
    ops.vids(i).exposureDur = exposureDur;
    ops.vids(i).Fs = FsPerColor;
end
%% for single color imaging
% ops.vids = ops.vids(1);
% ops.vids.frameMod = [1 0];


% create ROI
if makeROI == 1
    theseFiles = generateFileList(ops, 1);
    
    imageForROI = mean(loadTiffStack(theseFiles{1}, 'tiffobj',0),3); %make purple and blue image separately??
    imagesc(imageForROI);
    axis equal tight;
    colormap(gray);
    caxis(prctile(imageForROI(:),[1 99]));
    %colormouse;%change the colormap range interactively using the mouse
    title('left click to put anchor points., double click to finish');
    roiAhand = images.roi.AssistedFreehand;
    draw(roiAhand);
    ops.roi = createMask(roiAhand);
    close;
    
    %22/7/20 for use later
    thisROI=ops.roi;
    save(fullfile(rawDataDir, mouseName,'thisROI.mat'),'thisROI');
elseif makeROI == 0 %reuse ROI created previously.
    load(fullfile(rawDataDir, mouseName,'thisROI.mat'));
    ops.roi = thisROI; 
elseif makeROI == -1
    ops.roi = [];
end


%where SVD data is saved
ops.localSavePath = fullfile(rootDrive, 'data', mouseName, thisDate);

if ~exist(ops.localSavePath,'dir')
    mkdir(ops.localSavePath); %7/5/20
end
cd(ops.localSavePath)
save ops.mat ops


%% run SVD, save results to server, erase original data
pipelineHereKT();


%% save this script
originalFileNameAndLocation=[mfilename('fullpath')  '.m'];
suffix = ['_' mouseName '_' thisDate '_' num2str(thisSeries) '_' num2str(expNums)];
newFileLocation = fullfile(rawDataDir, mouseName, [thisDate '_' num2str(thisSeries)]);
[~, filename] = fileparts(originalFileNameAndLocation);
newFileNameAndLocation=fullfile(newFileLocation, [filename suffix '.m']);
A = exist(newFileNameAndLocation,'file');
if (A~=0)
    warning('Backup already exists for the current version')
else
    copyfile(originalFileNameAndLocation, newFileNameAndLocation);
end
