function [mDelay, mDelaySec] = phase2Delay(ANGLEMAPS, stimDur)
%[mDelay, mDelaySec] = phase2Delay(ANGLEMAPS, stimDur)

negative = ANGLEMAPS<0; %10/12
ANGLEMAPS(negative) = ANGLEMAPS(negative)+2*pi; %10/12

mDelay = 0.5*(ANGLEMAPS(:,:,1)+ANGLEMAPS(:,:,2)); %[rad]
%mDelay = circ_mean(ANGLEMAPS,[],3);
mDelay(mDelay<0) = mDelay(mDelay<0)+pi; %[0 2pi]

if nargin==2
    mDelaySec = stimDur/2/pi*mDelay; %[s]
else
    mDelaySec = [];
end