function [eyeFrameTimes, timeDiff] = grabEyeFrameTimes2P( expt, eyeRawPath, TSRawPath)
% [ frameTimes, timeDiff ] = grabEyeFrameTimes( expt, expr, TSDir)
% returns time of eye cam frame times in ThorSync time
% 25/6/20 created from MouseEyeTrack\+et\getFrameTimes

thisDate = expt.expDate(1:10);
thisSeries = str2num(expt.expDate(12:end));
DateSeriesNum = str2num([thisDate([1:4 6:7 9:10]) num2str(thisSeries)]);

TSFieldName = 'camStrobe';

%thr = 0.35; % Set this lower as in some cases there are slow changes in phd amplitude. Check this if the number of frames seems incorrect
thr = 0.5;

if nargin < 3
    TSRawPath = fullfile('D:\thorimagedata\', expt.subject, expt.expDate, num2str(expt.expNum));
    % will be replaced with:
    %     TSRawPath = fileparts(dat.expFilePath(expt.subject, thisDate, thisSeries, ...
    %         expt.expNum, '2p_raw','master'));
end
if nargin < 2
    %     eyeRawDir = 'D:\bhvCamData\';
    eyeRawPath = fileparts(dat.expFilePath(expt.subject, thisDate, thisSeries, ...
        expt.expNum, 'eyetracking_raw','master'));
end


%% retrieve exposure times from ThorSync
%fnam = fullfile(TSDir, subject, expDate, num2str(expr), 'Episode001.h5');
fnam = fullfile(TSRawPath, 'Episode001.h5');
[syncDataOut] = LoadSyncEpisodeFunction(fnam);

switch TSFieldName
    case 'camStrobe'
        tgtData = syncDataOut.camStrobe;
        
        above = tgtData>thr;
        deltas = [0; diff(above)];
        
        idx = find(deltas==1); %onset of the signal going HIGH
        
        nTotalFrames = length(idx);
end


frames = 1:nTotalFrames;

tsTimes = syncDataOut.time(idx(frames));
fprintf('There are %d frames in the ThorSync file\n', length(tsTimes));
% plot(syncDataOut.time, tgtData, times, 0, 'o');

if length(tsTimes) < 10 %4/12/20
    disp('camStrobe signal not saved properly. Use CaptureActive instead.');
    tgtData = syncDataOut.CaptureActive;
    
    above = tgtData>thr;
    deltas = [0; diff(above)];
    
    idx = find(deltas==1); %onset of the signal going HIGH
    
    tsTimes = syncDataOut.time(idx(1)) + 0.05;
end

%% retrieve frame times recorded by the camera & et.listen
%eyeRawPath = fullfile(eyeRawDir, expt.subject, num2str(DateSeriesNum), num2str(expt.expNum));
fileprefix = [num2str(DateSeriesNum) '_'  num2str(expt.expNum) '_'  expt.subject '_eye']; %or .mj2
load(fullfile(eyeRawPath,[fileprefix '.mat'])); %eyeLog
vr = VideoReader(fullfile(eyeRawPath, [fileprefix '.avi']));
videoFrames = vr.NumberOfFrames;

fprintf('There are %d frames in the video file\n', videoFrames);
fprintf('There are %d timestamps in the log file\n', length(eyeLog.TriggerData));

eyeFrameTimes = nan(videoFrames, 1);
for iFrame = 1:videoFrames
    if iFrame<=length(eyeLog.TriggerData)
        eyeFrameTimes(iFrame) = datenum(eyeLog.TriggerData(iFrame).AbsTime)*(24*60*60);
    else
        eyeFrameTimes(iFrame)=NaN;
    end
end

%idx = 1:min(length(tsTimes),length(eyeLog.TriggerData));
%timeDiff = median(eyeFrameTimes(idx)) - median(tsTimes(idx));
timeDiff = eyeFrameTimes(1) - tsTimes(1);


eyeFrameTimes = eyeFrameTimes - timeDiff;

