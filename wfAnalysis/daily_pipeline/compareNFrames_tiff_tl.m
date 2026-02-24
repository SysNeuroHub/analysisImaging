%  clear all
% addpath 'C:\Documents\git\analysisImaging'
% setPath_analysisImaging;
% 
% cd('C:\Documents\git\analysisImaging\wfAnalysis\daily_pipeline');
% load('amberRedOps.mat');
% % For blue/purple alternate
% %load('bluePurpleOps.mat'); 
% % For purple only
% % load('C:\Users\Experiment\Documents\MATLAB\purpleOps.mat')
% mouseName = 'Leekuanyew';%'susanoo';
% thisDate = '2026-02-19';%'2024-11-22'; %[datestr(now,'yyyy-mm-dd')];  
% thisSeries = 1;
% expNums = [1:2];
% thisDateSeries = [thisDate, '_' num2str(thisSeries)];
% 
% rootDrive = 'D:\svdinput';
% rawDataDir = 'D:\Subjects'; %local temporary storage
% 
% ops.mouseName = mouseName;
% ops.thisDate = thisDateSeries;
% ops.fileBase = fullfile(rawDataDir, mouseName, thisDateSeries);
% ops.statusDestination = fopen('test.txt','w');
% ops.rigName = 'alloptrig';%'wfrig'; %used in determineTimelineAlignments

function [nFr_tif, nFr_tl, status] = compareNFrames_tiff_tl(ops, expNums)
% ops.mouseName = mouseName;
% ops.thisDate = thisDateSeries;
% ops.fileBase = fullfile(rawDataDir, mouseName, thisDateSeries);
% ops.statusDestination = fopen('test.txt','w');
% ops.rigName = 'alloptrig';%'wfrig'; %used in determineTimelineAlignments

thisDate = ops.thisDate(1:10); %8/5/20
thisSeries = ops.thisDate(12:end); %8/5/20
for e = 1:length(expNums)
    expRefs{e} = dat.constructExpRef(ops.mouseName, thisDate, thisSeries, expNums(e));
end


nFr_tif = [];
nFr_tl = [];
for e = 1:length(expRefs)
%% tiff
    ops.inclExpList = expNums(e);
    theseFiles = generateFileList(ops);
    [nFr_tif(e), nFrPerFile] = getNFramesFromTifFiles(theseFiles, ops.statusDestination);


%% timeline    
    timelinePath = dat.expFilePath(expRefs{e}, 'timeline', 'master');
    if exist(timelinePath)
        load(timelinePath)
        strobeTimes = getStrobeTimes(Timeline, ops.rigName);
        theseStrobeNumbers = 1:numel(strobeTimes);
        inclStrobes = 1:numel(theseStrobeNumbers);
        % inclStrobes = mod(theseStrobeNumbers, ops.frameMod(1))==ops.frameMod(2);
        nFr_tl(e) = numel(strobeTimes(inclStrobes));
        % allT{e} = strobeTimes(inclStrobes);
    else
        disp(['NOT found ' timelinePath]);
    end

    status(e) = (nFr_tif(e)==nFr_tl(e));
    
    disp(['Exp' num2str(expNums(e)) ', tif: ' num2str(nFr_tif(e)) ', timeline: ' num2str(nFr_tl(e))]);
end
