function [ times] = grabExposeTimes( expt, expr, frames, TSDir, TSFieldName)
% [ times ] = grabExposeTimes( expt, expr, frames, TSDir)
% returns times of TTL signals going from LOW to HIGH
% for 2p frame time (default), return times of each plane as a cell (planes x 1)
% for other TSFieldNames ,return an array
%10/2/20 created from grabFrameTimes for thorsync version

subject = expt.subject;
expDate = expt.expDate;



if nargin < 4
    TSDir = 'D:\thorimagedata\'; %where ThorSync data is saved or copied
end

if nargin < 3
    frames = 'all';
end

%timelineFile = strcat(expDate,'_',num2str(expr),'_',subject,'_','Timeline.mat');
fnam = fullfile(TSDir, subject, expDate, num2str(expr), 'Episode001.h5');
[syncDataOut] = LoadSyncEpisodeFunction(fnam);
%OR [FTOSignalValid, stimulusSignalValid]=filterFTObyStimulus(fnam) ??

%ind = strcmp({Timeline.hw.inputs.name}, 'neuralFrames');

switch TSFieldName
    case 'slExposeClock' %5/2/20 ... will be separated into some different function
        tgtData = syncDataOut.slExposeClock; %or FTOSignalValid?
        
        phd = tgtData; % Helps with noisy photodiode signal
        
        thr = 0.35; % Set this lower as in some cases there are slow changes in phd amplitude. Check this if the number of frames seems incorrect
        
        above = phd>thr;
        deltas = [0; diff(above)];
        
        idx = find(deltas==1);
        
        nTotalFrames = length(idx);
end

if isequal(frames, 'all')
    frames = 1:nTotalFrames;
else
    frames = frames(frames<=nTotalFrames);
end

times = syncDataOut.time(idx(frames));






