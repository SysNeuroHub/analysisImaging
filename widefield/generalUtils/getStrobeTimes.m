function [strobeOnTimes, strobeOffTimes, strobeDurs] = getStrobeTimes(Timeline, rigName)
% [strobeOnTimes, strobeOffTimes, strobeDurs] = getStrobeTimes(Timeline, rigName)
% Returns the times in Timeline coordinates of every camera exposure 
%
% Inputs:
%   Timeline: timeline data
%   rigName: is a special parameter that identifies which Timeline data object should be used for this, and what threshold to apply to detect events. See the switch/case block below. 
%
% Outputs:
%   strobeOnTimes: times when strobe turns on
%   strobeOffTimes: times when strobe turns off
%   strobeDurs: durations when strobe is  on

switch rigName
    case 'bigrig'
        strobeName = 'cam2';
        strobeThresh = 2;
    case 'bigrig2'
        strobeName = 'cam2';
        strobeThresh = 2;    
    case 'bigrig1'
        strobeName = 'cam1';
        strobeThresh = 2;
    case 'kilotrode'
        strobeName = 'pcoExposure';
        strobeThresh = 2;    
    case 'wfrig' %8/5/20
        strobeName = 'camExposure';
        strobeThresh = 2;
    otherwise
        error('getStrobeTimes doesn''t recognize rig name %s', rigName);
end

strobesNum = find(strcmp({Timeline.hw.inputs.name}, strobeName));

ts = Timeline.rawDAQTimestamps;
strobes = Timeline.rawDAQData(:,strobesNum);

strobeOnSamps = find(strobes(1:end-1)<strobeThresh & strobes(2:end)>=strobeThresh);
strobeOffSamps = find(strobes(1:end-1)>strobeThresh & strobes(2:end)<=strobeThresh);

strobeOnTimes = ts(strobeOnSamps);
strobeOffTimes = ts(strobeOffSamps);

%% strobe duration 15/5/20
nSamps = min(length(strobeOnSamps), length(strobeOffSamps));
strobeDurs = strobeOffTimes(1:nSamps) - strobeOnTimes(1:nSamps);

