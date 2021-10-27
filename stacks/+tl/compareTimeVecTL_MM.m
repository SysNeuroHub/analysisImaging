function timeVecTL = compareTimeVecTL_MM(timeVecTL, timeVecMM, timeVecStTL, timeVecEnTL, TTLStTL, timeInfo)

%% compare memorymap.timevec (after fixing missed frames) and timeline.timevec
%for pco.edge camera

% 2015-9-17 DS fixed condition for length of timeInfo.elim

%% 2014-08-06 DS created from buildSingleTrStack.m

if (length(timeInfo.elim)>=1) && ((timeInfo.elim(1) == 1) || (timeVecStTL(1) < TTLStTL(1))) %include '=' ?. 23/4/14 DS
    timeVecTL = timeVecTL(2:end);
    disp('First frame eliminated (TimeLine)');
end

if (length(timeInfo.elim)>=2) && (timeInfo.elim(2) == 1)
    %do nothing to timeVecTL
end

while length(timeVecTL) > length(timeVecMM)
    warning('timeVec length does not match between timeline and stackset. Removing 1st frame (TimeLine)..');
    
    timeVecTL = timeVecTL(2:end);%just heuristics
end

if length(timeVecTL) > length(timeVecMM)
    error('length(timevec of timeline) > length(timevec of mmapfile)');
elseif length(timeVecTL) < length(timeVecMM)
    error('length(timevec of timeline) < length(timevec of mmapfile)');
end
