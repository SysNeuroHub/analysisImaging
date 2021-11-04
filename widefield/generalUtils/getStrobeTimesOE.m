function [strobeOnTimes, strobeOffTimes, strobeDurs] = getStrobeTimesOE(jsonFile, rigName)
% [strobeOnTimes, strobeOffTimes, strobeDurs] = getStrobeTimes(Timeline, rigName)
% Returns the times in Timeline coordinates of every camera exposure 
%
% Inputs:
%   jsonfile: fullpath to jsonfile OpenEphys data
%   rigName: is a special parameter that identifies which OpenEphys data object should be used for this, and what threshold to apply to detect events. See the switch/case block below. 
%
% Outputs:
%   strobeOnTimes: times when strobe turns on
%   strobeOffTimes: times when strobe turns off (NOT YET IMPLEMENTED)
%   strobeDurs: durations when strobe is  on (NOT YET IMPLEMENTED)

switch rigName
    case 'acuteMarmo'
        camStrobeCh = 7;
    otherwise
        error('getStrobeTimesOE doesn''t recognize rig name %s', rigName);
end

if ~exist(jsonFile, 'file')
    error(['DONT EXIST ' jsonFile]);
end

C = load_open_ephys_binary(jsonFile, 'continuous', 1);
%TODO: retrieve photodiode signal
fraxis = double(C.Timestamps); %this is frame number not time

%retrieve trial onset/offset
E = load_open_ephys_binary(jsonFile, 'events', 1);
srate = E.Header.sample_rate;
fr0 = fraxis(1);
taxis=1/srate*(fraxis-fr0); %[s]

theseEvents = find(E.ChannelIndex == camStrobeCh);
cam_ev = E.Data(theseEvents);
cam_fr = E.Timestamps(theseEvents)-fr0;

strobeOnTimes = taxis(cam_fr);
strobeOffTimes = []; %TO BE IMPLEMENTED

%% strobe duration 15/5/20
%nSamps = min(length(strobeOnSamps), length(strobeOffSamps));
%strobeDurs = strobeOffTimes(1:nSamps) - strobeOnTimes(1:nSamps);
strobeDurs = []; %TO BE IMPLEMENTED

