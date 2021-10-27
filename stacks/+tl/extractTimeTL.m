function [stIdx, enIdx] = extractTimeTL(Timeline, animal, iSeries, iExp, iStim, iTr)
%[stIdx, enIdx] = extractTimeTL(Timeline, animal, iSeries, iExp, iStim, iTr)
%extracts specified trial start/end time-indexes in timeline data of each trial

% 28-04-14 DS created
% 05-06-14 DS modified. Use TTL pulse from VS, when photoDiode is not recorded
% 31-05-15 DS iSeries can be string rather than number
% 24-07-15 DS separated out visStimSt (3rd output)
% 29-06-16 DS iSeries can be string
% 15-01-17 LFR if ~ischar(iSeries), convert iSeries to string


% TO DO: select photodiode or photosensor automatically using screenInfo ??


if ~ischar(iSeries)
    iSeries = num2str(iSeries);
end

margin = round(0.5*Timeline.hw.daqSampleRate); %10/8/2014 modefied
%margin size of time steps to allow data acquisition before and after UDP message
%find time when specified trial is started/ended
stString = sprintf('StimStart %s %s %d %d %d', animal, iSeries, iExp, iTr, iStim);
for ee = 1:length(Timeline.mpepUDPEvents)
    if any(strfind(Timeline.mpepUDPEvents{ee}, stString))
        stEvent = ee;
        stTime = Timeline.mpepUDPTimes(stEvent);
        break;
    end
end
if ~exist('stTime','var')
    error(['tl.extractTimeTL: ' stString ' NOT FOUND in Timeline']);
end

enString = sprintf('StimEnd %s %s %d %d %d', animal, iSeries, iExp, iTr, iStim);
for ee = 1:length(Timeline.mpepUDPEvents)
    if any(strfind(Timeline.mpepUDPEvents{ee}, enString))
        enEvent = ee;
        enTime = Timeline.mpepUDPTimes(enEvent);
        break;
    end
end
if ~exist('enTime','var')
    error(['tl.extractTimeTL: ' enString ' NOT FOUND in Timeline']);
end

%hack for M150219_SD
if strcmp(animal, 'M150219_SD') && enTime-stTime < 3;
    stTime = stTime - 10;
    disp('hack in extractTimeTL for M150219_SD');
end

[~, stIdx] = min(abs(Timeline.rawDAQTimestamps - stTime));
[~, enIdx] = min(abs(Timeline.rawDAQTimestamps - enTime));

stIdx = stIdx - margin;
enIdx = enIdx + margin;

if stIdx < 1 stIdx = 1; end
if enIdx > length(Timeline.rawDAQTimestamps) enIdx = length(Timeline.rawDAQTimestamps); end

