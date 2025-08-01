
clear all
setPath_analysisImaging;

%cd('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\master\Analysis\wfAnalysis\daily_pipeline');
cd('C:\Users\dshi0006\git\analysisImaging\wfAnalysis\daily_pipeline');
load('amberRedOps.mat');
mouseName = 'artemis';
thisDate = '2025-05-19';%'2024-11-22'; %[datestr(now,'yyyy-mm-dd')];
thisSeries = 1;
expNums = 1;%[1:4];
hwbinning = 1; %automatically retrieve this from thorcam header??
magnification = .5;
makeROI = 0; %1: make ROI and save thisROI.mat, 0: use already saved ROI from the save subject (thisROI.mat), -1: use all pixels
doRegistration = 0;%1; %15/10/20



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

%% load first record
theseFiles = generateFileList(ops, 1);

theseImages = loadTiffStack(theseFiles{1}, 'tiffobj',0); %make purple and blue image separately??
firstAmber = theseImages(:,:,1);


%% load ROI
try
    load(fullfile(rawDataDir, mouseName,'thisROI.mat'),'thisROI');
catch err
    imageForROI = mean(theseImages,3);
    imagesc(imageForROI);
    axis equal tight;
    colormap(gray);
    caxis(prctile(imageForROI(:),[1 99]));
    %colormouse;%change the colormap range interactively using the mouse
    title('left click to put anchor points., double click to finish');
    roiAhand = images.roi.AssistedFreehand;
    draw(roiAhand);
    thisROI = createMask(roiAhand);
    close all;
end



imagesc(firstAmber.*uint16(thisROI));

meanIntensity = mean(firstAmber(thisROI));

colormap(gray);
axis equal tight off;
colorbar;
caxis([0 500])
title([mouseName ' ' thisDate ', mean intensity: ' num2str(meanIntensity)])
screen2png([mouseName '-' thisDate]);
close all

