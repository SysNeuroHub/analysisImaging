function map_vf = phase2vfMap(ANGLEMAPS, vfRange, blankDur, stimDur, delay)
%map_vf = phase2vfMap(ANGLEMAPS, vfRange, blankDur, stimDur, delay)
%assumption: ANGLEMAPS(:,:,2) is response to increasing degree (rightward or upward)

if nargin < 5
     mDelay = 0.5*(ANGLEMAPS(:,:,1)+ANGLEMAPS(:,:,2)); %[rad]
    mDelay(mDelay<0) = mDelay(mDelay<0)+pi; %[0 2pi]
    
    delay = mean(mDelay(:)); %[0 2pi]
end

ANGLEMAPS = ANGLEMAPS - delay;
ANGLEMAPS(ANGLEMAPS < 0) = ANGLEMAPS(ANGLEMAPS < 0) + 2*pi; %[0 2pi]

map_phase = (ANGLEMAPS(:,:,2) - ANGLEMAPS(:,:,1));

doublePhaseRange = [-2*pi 2*pi*stimDur/(blankDur+stimDur)];%in theory this should be used but
%doublePhaseRange = [-pi pi*stimDur/(blankDur+stimDur)];%in practice this must be correct
map_vf = diff(vfRange)/diff(doublePhaseRange)*(map_phase - doublePhaseRange(1)) + vfRange(1); %[deg]
