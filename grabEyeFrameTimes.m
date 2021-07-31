function [times] = grabEyeFrameTimes( expt, expr, frames, rig)
% [times] = grabEyeFrameTimes( expt, expr, frames, rig)

if strcmp(rig, '2p')
    [times] = grabEyeFrameTimes2P( expt, expr, frames, TSDir);
elseif strcmp(rig, 'wf')
    [times] = grabEyeFrameTimesWF( expt, expr, frames, TSDir);
end    