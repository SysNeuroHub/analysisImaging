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
setPath_analysisImaging;

%cd('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\master\Analysis\wfAnalysis\daily_pipeline');
cd('C:\Users\dshi0006\git\analysisImaging\wfAnalysis\daily_pipeline');
load('amberRedOps.mat');
% For blue/purple alternate
%load('bluePurpleOps.mat'); 
% For purple only
% load('C:\Users\Experiment\Documents\MATLAB\purpleOps.mat')
mouseName = 'test';%'susanoo';
thisDate = '2025-05-08';%'2024-11-22'; %[datestr(now,'yyyy-mm-dd')];  
thisSeries = 3;
expNums = 1;%[1:4];
hwbinning = 1; %automatically retrieve this from thorcam header??
magnification = .5; 
makeROI = -1; %1: make ROI and save thisROI.mat, 0: use already saved ROI from the save subject (thisROI.mat), -1: use all pixels
doRegistration = 0;%1; %15/10/20


%where vidXraw.dat and vidXreg.dat are created (subsequently moved to the data server)
rootDrive = 'C:\svdinput'; %NG ... too small
%rootDrive = 'E:\svdinput';
%will be used in;
%ops.localSavePath

%where raw data is temporally downloaded must be under
%rawDataDir/(animal)/(session)/(expNum)
%rawDataDir = '\\vault-v2.erc.monash.edu.au\MNHS-dshi0006\Subjects';%if the raw data is already uploaded to the server
rawDataDir = 'M:\Subjects'; %market server ... too slow to load
%rawDataDir = 'E:\Subjects'; %local temporary storage


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

% retrieve camera sampling rate and exposure duration from Timeline
timelinePath = dat.expFilePath(expRefs{1}, 'timeline', 'master');
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




%%
% mouseName = 'Erlanger';
% thisDate = datestr(now, 'yyyy-mm-dd');
% expNums = [2 3];
% 
% 
% addpath(genpath('C:\Users\Experiment\Documents\GitHub\widefield'))
% addpath(genpath('\\zserver.cortexlab.net\Code\Rigging\main'));
% addpath(genpath('\\zserver.cortexlab.net\Code\Rigging\cb-tools'));
% 
% s = svdVid.listen;
% 
% s.ops.mouseName = mouseName;
% s.ops.thisDate = thisDate;
% 
% s.wizard
% 
% s.ops.vids(1).fileBase = fullfile('J:\data', mouseName, thisDate);
% s.ops.vids(2).fileBase = fullfile('J:\data', mouseName, thisDate);
% s.ops.vids(1).rigName = 'kilotrode';
% s.ops.vids(2).rigName = 'kilotrode';
% s.ops.useGPU = true;
% 
% for e = 1:length(expNums)
%     clear ed
%     ed.expRef = dat.constructExpRef(mouseName, thisDate, expNums(e));
%     s.addExpDat(ed)
% end
% 
% s.ops.localSavePath = fullfile('J:\data', mouseName, thisDate);
% 
% cd(fullfile('J:\data', mouseName, thisDate))
% ops = s.ops;
% save ops.mat ops
% pipelineHereKT