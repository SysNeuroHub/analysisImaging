function vsStTime = visStimStart(Timeline, animal, iSeries, iExp, iStim, iTr, ...
    photodiodeTh)
% vsStTime = visStimStart(Timeline, animal, iSeries, iExp, iStim, iTr) 
% returns stimulus onset time from the onset of recording of Timeline
%
% vsStTime = visStimStart(Timeline, animal, iSeries, iExp, iStim, iTr, ...
%    photodiodeTh) allows to set threshold for photodiode (default:1)

if nargin < 7
    photodiodeTh = 1; % 10/8/2014
end

%find time of visual stimulus start
[stIdx, enIdx] = tl.extractTimeTL(Timeline, animal, iSeries, iExp, iStim, iTr);

[dataSequence] = tl.extractSequenceTL(Timeline, stIdx, enIdx, 'photoDiode');

if max(dataSequence) > photodiodeTh
    
    vsStIdx = min(find(diff(dataSequence) > photodiodeTh)) + stIdx + 1;%8/8/14 added +1
    if isempty(vsStIdx)
        warning('Photodiode timing was not detected. Trying broader detection..');
        
        %if rising is slow...should be done more neatly 19/5/14 DS
        vsStIdx = min(find(dataSequence > mean(dataSequence))) + stIdx;
        
        %if photodiode signal is not acquired, vsStIdx ~= stIdx
        
        if isempty(vsStIdx)
            error('Photodiode timing was not detected!');
        end
        
    end
else
    %added DS on 5/6/14
    warning('Photodiode signal was not properly recorded. Calcultate time from vs-TTL onset.');
    [dataSequence] = tl.extractSequenceTL(Timeline, stIdx, enIdx, 'syncEcho');
    aiTh = 0.5*max(dataSequence);
    vsStIdx = min(find(diff(dataSequence) > aiTh)) + stIdx;
end

vsStTime = Timeline.rawDAQTimestamps(vsStIdx);
