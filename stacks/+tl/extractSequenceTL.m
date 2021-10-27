function [dataSequence] = extractSequenceTL(Timeline, stIdx, enIdx, tgtString)
%[dataSequence] = extractSequenceTL(Timeline, stIdx, enIdx, tgtString)
% extracts data sequence from Timeline structure from StIdx to enIdx

tgtCh = find(strcmp({Timeline.hw.inputs.name}, tgtString));

if ~isempty(tgtCh)
    dataSequence = Timeline.rawDAQData(stIdx:enIdx, tgtCh);
else
    error(['''' tgtString ''' was not found in Timeline!']);
end