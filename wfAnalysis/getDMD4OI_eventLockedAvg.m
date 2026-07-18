function [DMD4OIsequence, tevent_DMD] = getDMD4OI_eventLockedAvg(image4OI, ...
laser, DMDidx, t, stimSequence, laserTimes, respWin)
%[DMD4OIsequence, tevent_DMD] = getDMD4OI_eventLockedAvg(image4OI, ...
% DMDOut_state, Timeline, t, stimSequence, stimTimes, nstim, respWin)
% returns event-triggered DMD image, with respect to respWin
% ... not happy with the inputs

% TODO: handling when respWin exceeds timeline time

nstim = numel(stimSequence.labels);

if size(laser,2) == 1; laser = laser'; end
[avgPeriEvent_lsr] = ...
    eventLockedAvg(laser, t, laserTimes.onset, stimSequence.seq, respWin);

if size(DMDidx,2) == 1; DMDidx = DMDidx'; end
[avgPeriEvent_DMD, tevent_DMD] = ...
    eventLockedAvg(DMDidx, t, laserTimes.onset, stimSequence.seq, respWin);

DMD4OIsequence = zeros(size(image4OI,1), size(image4OI,2),numel(tevent_DMD), nstim);
for istim = 1:nstim
    DMD4OIsequence(:,:,:,istim) = image4OI(:,:,round(avgPeriEvent_DMD(istim,1,:))) .* avgPeriEvent_lsr(istim,1,:);
end

