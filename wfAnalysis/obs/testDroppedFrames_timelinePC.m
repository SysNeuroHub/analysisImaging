addpath(genpath('C:\Users\Experiment\Documents\MATLAB\visbox'));
addpath(genpath('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\widefield_ns'));

%% load Timeline
load('C:\Users\Experiment\Documents\MATLAB\Data\Subjects\dummy_wf\2020-05-08_9\1\2020-05-08_9_1_dummy_wf_Timeline.mat');
tlExposeClockIdx = find(strcmp({Timeline.hw.inputs.name}, 'tlExposeClock'));
ExposeClockTimes = find(Timeline.rawDAQData(1:end-1,tlExposeClockIdx) <=2 & ...
    Timeline.rawDAQData(2:end,tlExposeClockIdx) > 2);
camExposureIdx = find(strcmp({Timeline.hw.inputs.name}, 'camExposure'));
camExposureTimes = find(Timeline.rawDAQData(1:end-1,camExposureIdx) <=2 & ...
    Timeline.rawDAQData(2:end,camExposureIdx) > 2);


%% load tiff files
cd('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\master\Analysis\wfAnalysis\daily_pipeline');
% For purple only
load('purpleOps.mat')
mouseName = 'dummy_wf';
thisDate = '2020-05-08';%[datestr(now,'yyyy-mm-dd')];
thisSeries = 9;
expNums = 1;
theseFiles = {'C:\ThorCamData\dummy_wf\2020-05-08_9\1.tif','C:\ThorCamData\dummy_wf\2020-05-08_9\1_1.tif'};



rootDrive = 'C:\svdinput';
rawDataDir = 'C:\ThorCamData'; 


%% modify ops for this experiment
thisDateSeries = [thisDate, '_' num2str(thisSeries)];
ops.mouseName = mouseName;
ops.thisDate = thisDateSeries;
%ops.inclExpList = '1';%'exp1_uncompressed';
ops.fileBase = rawDataDir;
ops.statusDestination = fopen('test.txt','w');
ops.userName = 'Daisuke';
ops = rmfield(ops,'emailAddress');
ops.rigName = 'wfrig'; %used in determineTimelineAlignments

for i = 1:length(ops.vids)
    ops.vids(i).fileBase = fullfile(rawDataDir, mouseName, thisDateSeries);
end

for e = 1:length(expNums)
    expRefs{e} = dat.constructExpRef(mouseName, thisDate, thisSeries, expNums(e));    
end
ops.expRefs = expRefs;

%where SVD data is saved
ops.localSavePath = fullfile(rootDrive, 'data', mouseName, thisDate);

mkdir(ops.localSavePath); %7/5/20
cd(ops.localSavePath)

addpath(genpath('C:\Users\Experiment\Documents\MATLAB\npy-matlab'));
addpath('C:\svdinput');

serverDir = '\\ad.monash.edu\home\User006\dshi0006\Documents\tempDataServer'; %7/5/20
%'\\lugaro.cortexlab.net\bigdrive\staging\';

%load ops.mat; % this must be present in the current directory
% diaryFilename = sprintf('svdLog_%s_%s.txt', ops.mouseName, ops.thisDate);
% diary(diaryFilename);
    
ops.localSavePath = pathForThisOS(ops.localSavePath);
for v = 1:length(ops.vids)
    ops.vids(v).fileBase = pathForThisOS(ops.vids(v).fileBase);
end
ops.fileBase = ops.vids.fileBase; %11/5/20

if ~exist(ops.localSavePath, 'dir')
    mkdir(ops.localSavePath);
end
save(fullfile(ops.localSavePath, 'ops.mat'), 'ops');    

%% load all movies into flat binary files, for each color
for v = 1:length(ops.vids)
    
    clear loadDatOps;
    
    ops.theseFiles = [];
    %theseFiles = generateFileList(ops, v);
    
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
    %check how frameNumbersFromStamp & frameNumbersWithinRec are used from
    %here
    %%sanity check
    % subplot(211);
    %     plot(dataSummary.frameRecIndex); hold on
    %     plot(dataSummary.frameFileIndex);
    %     legend('frameRecIndx','frameFileIndex');
    %     grid on;
    % subplot(212);
    % plot(dataSummary.imageMeans);
    
    fn = fieldnames(dataSummary);
    results(v).name = ops.vids(v).name;
    for f = 1:length(fn)
        results(v).(fn{f}) = dataSummary.(fn{f});
    end
    
    %save(fullfile(ops.localSavePath, 'results.mat'), 'results');
end

disp(['TL exposeclock:' num2str(length(ExposeClockTimes))]);
disp(['TL camExposure:' num2str(length(camExposureTimes))]);
disp(['TIFF frames:' num2str(length(dataSummary.imageMeans))]);

