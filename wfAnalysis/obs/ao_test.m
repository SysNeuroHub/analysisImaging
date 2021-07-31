%% this is a script to do the following:
% detect onset of each stimulation in thorsync time using photodiode signal
% detect wf frame time in thorsync time
% obtain stimulus triggered time courses
% compute Direction selectivity index of each cell

%TODO
% analysis for ORsequence

addpath(genpath('C:\npy-matlab'));
addpath(genpath('C:\Users\dshi0006\npy-matlab'));
addpath(genpath('C:\Users\Analysis\npy-matlab'));
addpath(genpath('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\master\Analysis\wfAnalysis\widefield'));
addpath(genpath('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\master\visbox\rigbox'));
addpath('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\master\Analysis');
addpath('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\master\Analysis\wfAnalysis');


%% experiment
%ppbox notation
expt.subject = 'dummy_wf';
expt.expDate = '2020-09-25_1';
expt.expNum = 8;
bklightCtrl = 1;

%% analysis
%highpassCutoff = 0.01; %[Hz]
%lowpassCutoff = []; %[Hz]
marginT = 1; %[s]
resizeS = 1; %spatial rescaling factor


resultSaveDir = fullfile('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\AnalysisResult\wf',...
    [expt.subject '_' expt.expDate '_' num2str(expt.expNum)]);
mkdir(resultSaveDir);


%% data location
thisDate = expt.expDate(1:10);
thisSeries = str2num(expt.expDate(12:end));
expPath = fileparts(dat.expFilePath(expt.subject, thisDate, thisSeries, expt.expNum, 'widefield','master'));
saveVpath = expPath;
mpepDir = dat.reposPath('main', 'master');

 
%% load mpep data
p = ProtocolLoad_wf(expt.subject,expt.expDate,expt.expNum); %3/6/20

stimSequence = getStimSequence(expt.subject, expt.expDate, expt.expNum);
%stimSequence.seq
%stimSequence.labels
%stimSequence.paramValues
if ~isfield(p,'pfiledurs')
    p.pfiledurs = p.pars(1,:)/10;
end

respWin = [-marginT min(p.pfiledurs)+marginT]; %31/3/20


%% detect onset of each stimulation in Timeline time using photodiode signal
getAllPhdFlips = 0;
expt = grabStimTimesWF(expt, getAllPhdFlips, [], [], bklightCtrl);


%% detect onset of analogue output
TLDir = fileparts(dat.expFilePath(expt.subject, thisDate, thisSeries, expt.expNum, 'timeline','master'));
filename_TL = sprintf('%s_%d_%s_Timeline.mat',expt.expDate, expt.expNum, expt.subject);
timeline_filename = fullfile(TLDir, filename_TL);
load(timeline_filename);

%taken from grabStimTimesWF.m photodiode part
 tltime = Timeline.rawDAQTimestamps';
   
ao_idx = strcmp({Timeline.hw.inputs.name}, 'ao0');
screenOnTh = mean(Timeline.rawDAQData(:,ao_idx)); %25/9/20
ao_trace_raw = Timeline.rawDAQData(:,ao_idx)';
ao_on = ao_trace_raw > screenOnTh;

ao_trace = zeros(size(ao_trace_raw));
ao_trace(ao_on)=1;

%this judges whether the ao signal flipped or continued
ao_flip = find((~ao_trace(1:end-1) & ao_trace(2:end)) | ...
    (ao_trace(1:end-1) & ~ao_trace(2:end)))+1; %time idx when signal flips 0 to 1 or 1 to 0 ??


ao = struct('timestamps',[],'values',[]);
ao.timestamps = tltime(ao_flip);
ao.values = ao_trace(ao_flip)';

ao_offsets = ao.timestamps(ao.values == 0); %maybe missing/adding last frame of each repeat
ao_onsets = ao.timestamps(ao.values == 1);


%% compare photodiode v analogue output onsets
tdiff = expt.stimTimes.onset - ao_onsets;
histogram(1e3*tdiff,10);
std(1e3*tdiff)
mean(1e3*tdiff)
