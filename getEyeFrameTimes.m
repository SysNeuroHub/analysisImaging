function [frameTimes, tl_flag] = getEyeFrameTimes(varargin)

% This function will return the frametimes of the eye-tracking movie
% The frame times will be aligned to the Timeline time axis
% If corresponding Timeline file does not exist, frame times will be
% relative only.
%
% different possible ways to call this function:
%
% frameTimes = getFrameTimes(animal, series, exp)
%              animal, series, experiment - as used in mpep 
%              (series and exp are numericals)
%              frameTimes - timestamps of frames in seconds, 
%              aligned to the Timeline time axis
%
% 25/6/26 DS created from MouseEyeTrack/+et/getFrameTimes
%
% TODO: handling when nFrames from movie, timeline or log do not match
%
%this code relies on time difference of the same UDP events recorded in eyeLog and Timeline

useVideoTimeStamp = 0; 
%if 1, use time stamps in video as frameTimes
%if 0, use camStrobe in Timeline as frameTimes

if nargin == 3
    animal = varargin{1};
    series = varargin{2};
    exp = varargin{3};
    if isnumeric(exp)
        exp = num2str(exp);
    end
    %str = sprintf('ExpStart %s %s %d', animal, series, exp);
    %info = dat.mpepMessageParse(str);
    info.expRef = dat.constructExpRef(animal, series(1:10), str2double(series(12:end)), str2double(exp));
%     eyeFolder = fileparts(dat.expFilePath(animal, series(1:10), series(12:end), num2str(exp), 'widefield','master'));
    fullNames = dat.expFilePath(info.expRef, 'eyetracking');
%     % the second cell is the server location
     [eyeFolder ,eyeFileStem, ~] = fileparts(fullNames{2});
%     eyeFolder = ['\\zserver' eyeFolder(11:end)];
    % loading the eyeLog data
    warning off % there is always an annoying warning about a videoinput object
    load(fullfile(eyeFolder, eyeFileStem),'eyeLog');  %eyeLog
    warning on

    fullNames = dat.expFilePath(info.expRef, 'Timeline','master');%20/6/18
    % the second cell is the server location
    % loading the Timeline data
    load(fullNames,'Timeline');

elseif nargin == 2
    % the two arguments are the two filenames - eye and TL
elseif nargin == 1
    % only eye-tracking file is supplied, will figure out alone where is
    % the Timeline file
elseif nargin == 0
    % will open uigetfile(), so that user will be able to choose the video
    % file
    startPath = '\\zserver2\Data\EyeCamera\';
    [Filename, eyeFolder] = uigetfile('*.*', 'Choose an eye-tracking data file', startPath);
    [~, eyeFileStem, ~] = fileparts(Filename);
    % loading the eye data
    warning off % there is always an annoying warning about a videoinput object
    load(fullfile(eyeFolder, eyeFileStem),'eyeLog');
    warning on
    
    % figuring out where the Timeline data sits
    und_idx = strfind(eyeFileStem, '_eye');
    expRef = eyeFileStem(1:und_idx(end)-1);
end
    
vReader = VideoReader(fullfile(eyeFolder, eyeLog.loggerInfo.Filename));
nFrames_vid = vReader.NumFrames; %NumberOfFrames
nFrames_log = length(eyeLog.TriggerData);
fprintf('There are %d frames in the video file\n', nFrames_vid);
fprintf('There are %d timestamps in the log file\n', nFrames_log);

% find the last ExpStart-ExpEnd pair
endInd = [];
startInd = [];
for iEvent = length(eyeLog.udpEvents):-1:1
    if strfind(eyeLog.udpEvents{iEvent}, 'ExpEnd')
        endInd = iEvent;
        break;
    end
end
for iEvent = (endInd-1):-1:1
    if strfind(eyeLog.udpEvents{iEvent}, 'ExpStart')
        startInd = iEvent;
        break;
    end
end

udpEvents = eyeLog.udpEvents(startInd:endInd);
udpEventTimes = eyeLog.udpEventTimes(startInd:endInd);

nEvents = length(udpEventTimes);

if exist('Timeline','var')
    tl_flag = 1;
else
    tl_flag = 0;
    warning('Timeline was not found!!!');
end

% if  tl_flag && ~isequal(nEvents, Timeline.mpepUDPCount)
%     warning('Number of UDP events logged by Timeline and by EyeCamera is different. Something must be wrong!');
% end

% converting absolute times to times in seconds
eyeTimes = nan(nEvents, 1);
for iEvent=1:nEvents
    eyeTimes(iEvent) = datenum(udpEventTimes{iEvent})*(24*60*60);
end
frameTimes_vid = nan(nFrames_vid, 1);
for iFrame = 1:nFrames_vid
    if iFrame<=length(eyeLog.TriggerData)
        frameTimes_vid(iFrame) = datenum(eyeLog.TriggerData(iFrame).AbsTime)*(24*60*60);
        %        frameTimes(iFrame) = datetime(eyeLog.TriggerData(iFrame).AbsTime)*(24*60*60);
    else
        frameTimes_vid(iFrame)=NaN;
    end
end

if tl_flag
    tlTimes = Timeline.mpepUDPTimes(1:nEvents);
    nEvents2Discard = 2; % the first few event timing is unreliable
    idx = nEvents2Discard+1:nEvents;
    timeDiff = median(eyeTimes(idx)) - median(tlTimes(idx));
    frameTimes_vid = frameTimes_vid - timeDiff;    
    
    %% timestamps recorded in timeline
    strobetime = Timeline.rawDAQTimestamps';
    camStrobe_idx = strcmp({Timeline.hw.inputs.name}, 'camStrobe');
    %         camStrobe_flicker = max(Timeline.rawDAQData(:,camStrobe_idx)) - ...
    %             min(Timeline.rawDAQData(:,camStrobe_idx)) > 2;
    camStrobe_th = max(Timeline.rawDAQData(:,camStrobe_idx))/2;
    camStrobe_on = Timeline.rawDAQData(:,camStrobe_idx) < camStrobe_th; %[0 1]
    camStrobeEv = trace2Event(camStrobe_on, strobetime);
    frameTimes_tl = mean(camStrobeEv,2);
    nFrames_tl = numel(frameTimes_tl);
    fprintf('There are %d frames in Timeline (camStrobe) \n', nFrames_tl);%23/6/26
end

if nFrames_vid == nFrames_tl+1
    frameTimes_vid = frameTimes_vid(2:end);
    nFrames_vid = nFrames_vid - 1;
    disp('the first frame in eye video removed');
end

if nFrames_vid ~= nFrames_tl
    disp('#frames in eye video does NOT match #strobes in Timeline')
    %DO SOMETHING HERE
end
if nFrames_tl ~= nFrames_log
    %DO SOMETHING HERE
end
if nFrames_vid ~= nFrames_log
    %seems like this never happens
end

if useVideoTimeStamp
    frameTimes = frameTimes_vid;
else
    frameTimes = frameTimes_tl;
end

return;

%% =============some plotting for debugging purposes=========
% eyeTimes = eyeTimes - timeDiff;
% 
% figure
% stem(eyeTimes, ones(nEvents, 1), 'b');
% hold on;
% stem(tlTimes, ones(nEvents, 1), 'r:');
% legend('eyeUDPs', 'tlUDPs');
% xlabel('time [seconds]');
% title('UDP messges Timing (aligned)');
% 
% figure
% dd = diff(eyeTimes-tlTimes);
% plot(dd(nEvents2Discard+1:end))
% title('UDP timing jitter (eyeCamera - Timeline)');
% xlabel('UDP message number');
% ylabel('Time difference [sec]');
% 
% figure
% hist(dd(nEvents2Discard+1:end), 20);
% title('Time jitter histogram');
% xlabel('Time difference [sec]');
% 
% 
