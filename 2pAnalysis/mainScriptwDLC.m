%% this is a script to do the following:
% detect onset of each stimulation in thorsync time using photodiode signal
% detect bhv camra frame time in thorsync time
% obtain stimulus triggered time courses
% compute Direction selectivity index of each cell

%% TODO
%rename DLC file


%% experiment 
%ppbox notation
expt.subject = 'Copernicus';
expt.expDate = '2020-01-17';
expt.expNum = 5;

%% data location
TSDir = 'D:\thorimagedata\'; % = ops0.RootStorage where ThorSync data is saved or copied
DLCDir = 'D:\OutputDLC\';
mpepDir = 'Y:\'; %= ops0.mpepRootStorage

%% analysis
%iplane = 1;
respWin = [0 2]; %[s]
bodyName = 'eye_l';
coordsName = 'x';


%% load mpep data
p = ProtocolLoad(expt.subject,expt.expDate,expt.expNum, 'donotload', mpepDir);
stimSequence = getStimSequence(expt.subject, expt.expDate, expt.expNum, mpepDir); 
%stimSequence.seq
%stimSequence.labels
%stimSequence.paramValues


%% load DLC result (deeplabcut.analyze_videos(config_path, [video], save_as_csv=True))
% s2pname = sprintf('%s/%s/%s/%d/F_%s_%s_plane%d_proc.mat', DLCDir, ...
%         expt.subject, expt.expDate, expt.expNum, expt.subject, expt.expDate, iplane);
dlcname = 'D:\bhvCamData\stube_2020-02-03-091404-0000DeepCut_resnet50_2prig_test3Feb1shuffle1_30000.csv';
csvtable = readtable(dlcname);
%t2 = readmatrix(dlcname);%will be removed in a future release
%t3 = csvread(dlcname); %not recommended

traces = csvtable{3:end, 2:end}; % time x nameIdx
traces = cellfun(@str2double,traces);
bodyCoords = csvtable{1:2, 2:end}; %[(bodyparts, coords), entry]
bodyCoords = string(bodyCoords);
bodyparts = unique(bodyCoords(1,:));
coords = unique(bodyCoords(2,:));

%% use only specified bodyparts and coords names
nameIdx = intersect(find(strcmp(bodyCoords(1,:), bodyName)), ...
    find(strcmp(bodyCoords(2,:), coordsName)));
traces = traces(:,nameIdx);
ncells = length(nameIdx);


%% detect onset of each stimulation in thorsync time using photodiode signal
[ expt ] = grabStimTimes( expt, false, TSDir );
% expt.stimTimes(expt.expNum).onset
% expt.stimTimes(expt.expNum).offset
stimTimes = expt.stimTimes.onset; %[#stim onset x 1] time in thorsync time

%check number of stimulation
if length(stimTimes) ~= p.nstim * p.nrepeats
    error(['Detected #stim onset:' num2str(stimOnTime) '/ total #stim' num2str(p.nstim * p.nrepeats)]);
end

%% detect camera frame time in thorsync time
[ frameTimes ] = grabExposeTimes( expt, expt.expNum, 'all', TSDir, 'slExposeClock'); %[1 x 2pframes] 6/2/20

infoInXml = loadThorlabsExperimentXml(fullfile(TSDir, expt.subject, expt.expDate, ...
    num2str(expt.expNum)));
movStride = infoInXml.numZSlices; 
frameTimes = frameTimes(iplane:movStride:end); %maybe WRONG when triggered by TTL??

%% match thorsync and thorimage frames
%[FTOSignalValid, stimulusSignalValid]=filterFTObyStimulus(fnam) ??
if length(frameTimes) ~= size(traces,1)
    warning(['TI frames:' num2str(size(traces,1)) ', TS frames:' num2str(length(frameTimes))]);
    frameTimes = frameTimes(1: size(traces,1)); %hack for now
end

    

%% stimulus triggered time courses
[avgPeriEventV, winSamps] = eventLockedAvg(traces', frameTimes, stimTimes, stimSequence.seq, respWin);
% avgPeriEventV: nEventTypes x nCells x nTimePoints
% winSamps: labels for the time axis, relative to the event times
eLabel = unique(stimSequence.seq); %1st dimension of avgPeriEventV
eLabel_paramValue = stimSequence.paramValues(eLabel);

stimIdx = eLabel(~isnan(eLabel_paramValue));
blankIdx  = eLabel(isnan(eLabel_paramValue));
stimValue = eLabel_paramValue(stimIdx);
avgPeriEventV = avgPeriEventV(stimIdx, :, :);% - avgPeriEventV(blankIdx, :, :); %subtract by blank condition


