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


%% TODO
%apply the SVD to multipe exps at once
%fix ops.expRefs
%oephys file name to be renamed at the time of recording?
%retrieve camera sampling rate from the record
%timestamps not recorded

%dsbox,oi-tools,analysis-tools,analysisImaging, npy-matlab,neurostim, marmolab-stimuli
if exist('C:\Users\dshi0006\git','dir')
    addpath(genpath('C:\Users\dshi0006\git'));
else
    addpath(genpath('C:\git'));    
end

addpath(genpath('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\master\visbox'));
%dat.expFilePath, dat.constructExpRef
%addpath('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\master\visbox\rigbox\cb-tools\burgbox');
%getOr

clear all



workDir = 'C:\Users\dshi0006\git\analysisImaging\wfAnalysis\daily_pipeline';
cd(workDir);
load('intrinsicImagingOps.mat'); 

%% about experiment
subject = 'rat1gou';
expDate = '20211026';
% thisDate = '2021-01-13'; %[datestr(now,'yyyy-mm-dd')];  
thisSeries = 1; %tentative
expNums = 2:3;%[1:5];%[1:4];
makeROI = true;

%% about imaging
hwbinning = 2; %automatically retrieve this from thorcam header??
magnification = 1; %0.5
doRegistration = 1; %15/10/20
FsPerColor = 10;%[Hz] retrieve this from OI data/header
exposureDur = 1/FsPerColor; %[s]

%% about openEphys recording camera strobes
OErecDir = '\\storage.erc.monash.edu\shares\R-MNHS-Syncitium\Shared\MarmosetData';

%for multiple experiments, append multiple filenames as cell
DirBase_before{1} = fullfile(OErecDir, subject,'oephys','experiment1/recording1');



%where vidXraw.dat and vidXreg.dat are created (subsequently moved to the data server)
rootDrive = 'C:\svdinput';
%will be used in;
%ops.localSavePath

%where raw data is saved must be under
%rawDataDir/(animal)/(session)/(expNum)
%rawDataDir = '\\vault-v2.erc.monash.edu.au\MNHS-dshi0006\Subjects';%if the raw data is already uploaded to the server
rawDataDir = '\\storage.erc.monash.edu\shares\R-MNHS-Syncitium\Shared\MarmosetData'; 

%where SVD and summary data is saved
%dat.expFilePath(ops.mouseName, thisDate, thisSeries, expNums, 'widefield','master')

%where the raw and processed data is saved
%dat.expFilePath(ops.mouseName, thisDate, thisSeries, expNums, 'widefield_raw','master')


%% modify ops for this experiment
%thisDateSeries = [thisDate, '_' num2str(thisSeries)];
if ~makeROI
    if exist(fullfile(rawDataDir, subject,'thisROI.mat'),'file')
        load(fullfile(rawDataDir, subject,'thisROI.mat'),'thisROI');
        ops.roi = thisROI;
    else
        ops.roi = [];
    end
end

ops.mouseName = subject;
ops.thisDate = expDate;
%if using inclExpList, specify ops.fileBase
ops.inclExpList = expNums;
%where original video data is stored, where generateFileList looks for data
ops.fileBase = fullfile(rawDataDir, subject, 'imaging',expDate);
ops.statusDestination = fopen('test.txt','w');
ops.userName = 'Daisuke';
ops.rigName = 'acuteMarmo'; %used in determineTimelineAlignments
ops.doRegistration = doRegistration;
ops.useGPU = 1; %used only for registration?
ops.objectiveType = num2str(magnification); %0.5 / 0.8
ops.pixelSizeUM = 5/str2num(ops.objectiveType)*hwbinning; %added 15/5/20
ops.hasASCIIstamp = 0; %18/5/20

%this is used where to save the SVD results in the server
for e = 1:length(expNums)
    expRefs{e} = dat.constructExpRef(subject, expDate, thisSeries, expNums(e));    
end
ops.expRefs = expRefs;

% % % retrieve camera sampling rate and exposure duration from Timeline
% % timelinePath = dat.expFilePath(expRefs{1}, 'timeline', 'master');
% % load(timelinePath);
% % [strobeOnTimes, ~, strobeDurs] = getStrobeTimes(Timeline, ops.rigName);
% % exposureDur = median(strobeDurs); %[s]
% % FsPerColor = 1/median(diff(strobeOnTimes))/length(ops.vids);
% % ops.Fs = FsPerColor;%16/10/20
% % clear Timeline

for i = 1:length(ops.vids)    
    ops.vids(i).fileBase = ops.fileBase;
    ops.vids(i).exposureDur = exposureDur;
    ops.vids(i).Fs = FsPerColor; 
    ops.vids(i).thisDatPath = fullfile(ops.localSavePath, ['vid' num2str(i) 'raw.dat']);
end

%% create ROI ... to be used at the stage of SVD processing
if makeROI
    ops = addROItoOPS(ops);
    
    %22/7/20 for use later
    thisROI=ops.roi;
    save(fullfile(rawDataDir, subject, 'imaging','thisROI.mat'),'thisROI');
else
    load(fullfile(rawDataDir, subject, 'imaging','thisROI.mat'));
    ops.roi = thisROI;
end


%% where SVD data is saved temporally
ops.localSavePath = fullfile(rootDrive, 'data', subject, expDate);

if ~exist(ops.localSavePath,'dir')
    mkdir(ops.localSavePath); %7/5/20
end
cd(ops.localSavePath)
save ops.mat ops


%% rename openEphys files
expt.subject = subject;
expt.expDate = expDate;
expt.expNum = expNums;
renameOEphysFiles(expt, DirBase_before);

%% run SVD, save results to server, erase original data
pipelineHereKTns();




%%
% mouseName = 'Erlanger';
% expDate = datestr(now, 'yyyy-mm-dd');
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
% s.ops.thisDate = expDate;
% 
% s.wizard
% 
% s.ops.vids(1).fileBase = fullfile('J:\data', mouseName, expDate);
% s.ops.vids(2).fileBase = fullfile('J:\data', mouseName, expDate);
% s.ops.vids(1).rigName = 'kilotrode';
% s.ops.vids(2).rigName = 'kilotrode';
% s.ops.useGPU = true;
% 
% for e = 1:length(expNums)
%     clear ed
%     ed.expRef = dat.constructExpRef(mouseName, expDate, expNums(e));
%     s.addExpDat(ed)
% end
% 
% s.ops.localSavePath = fullfile('J:\data', mouseName, expDate);
% 
% cd(fullfile('J:\data', mouseName, expDate))
% ops = s.ops;
% save ops.mat ops
% pipelineHereKT