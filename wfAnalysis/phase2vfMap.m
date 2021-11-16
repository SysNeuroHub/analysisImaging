function map_vf = phase2vfMap(ANGLEMAPS, vfRange, blankDur, stimDur)
%map_vf = phase2vfMap(ANGLEMAPS, vfRange, blankDur, stimDur)

map_phase = (ANGLEMAPS(:,:,2) - ANGLEMAPS(:,:,1));

doublePhaseRange = [-2*pi 2*pi*stimDur/(blankDur+stimDur)];%in theory this should be used but
%doublePhaseRange = [-pi pi*stimDur/(blankDur+stimDur)];%in practice this must be correct
map_vf = diff(vfRange)/diff(doublePhaseRange)*(map_phase - doublePhaseRange(1)) + vfRange(1); %[deg]
