function [DMDIn_state, DMDOut_state] = getDMDState(Timeline, nPtn)
%  [DMDIn_state, DMDOut_state] = getDMDState(Timeline, nPtn)
% returns DMD pattern index at each timeline time
% nPtn: total number of patterns stored in .seq file
%     nPtn = parameters.Protocol.pars(2);

DMDInch = find(strcmp({Timeline.hw.inputs.name}, 'DMDIn'));
DMDOutch = find(strcmp({Timeline.hw.inputs.name}, 'DMDOut'));

DMDIn_raw = Timeline.rawDAQData(:,DMDInch);
DMDIn = DMDIn_raw > 2.5;
DMDOut_raw = Timeline.rawDAQData(:,DMDOutch);
DMDOut = DMDOut_raw > 2.5;

DMDIn_state = [0; mod_nonzero(cumsum(diff(DMDIn)>0), nPtn)];
DMDOut_state = [0; mod_nonzero(cumsum(diff(DMDOut)>0), nPtn)];
