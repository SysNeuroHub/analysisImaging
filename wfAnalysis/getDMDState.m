function [DMDIn_state, DMDOut_state] = getDMDState(Timeline, nPtn)
%  [DMDIn_state, DMDOut_state] = getDMDState(Timeline, nPtn)
% returns DMD pattern index at each timeline time
% nPtn: total number of patterns stored in .seq file
%     nPtn = parameters.Protocol.pars(2);

Vth = 2; %2.5 [V]
medfiltLength = 500; %50ms when timeline is sampled at 10kHz

DMDInch = find(strcmp({Timeline.hw.inputs.name}, 'DMDIn'));
DMDOutch = find(strcmp({Timeline.hw.inputs.name}, 'DMDOut'));

DMDIn_raw = Timeline.rawDAQData(:,DMDInch);
DMDIn = DMDIn_raw > Vth;
DMDOut_raw = Timeline.rawDAQData(:,DMDOutch);
DMDOut = DMDOut_raw > Vth;

DMDIn_state = medfilt1([0; mod_nonzero(cumsum(diff(DMDIn)>0), nPtn)], medfiltLength);
DMDOut_state = medfilt1([0; mod_nonzero(cumsum(diff(DMDOut)>0), nPtn)], medfiltLength);
end

function r = mod_nonzero(a,b)
%when a == b, return b instead of 0
r = mod(a-1, b) + 1;
end