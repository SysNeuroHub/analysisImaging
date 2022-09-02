function [taxis, fraxis] = getTAxisOE(jsonFile)
% taxis = getTAxisOE(jsonFile)

C = load_open_ephys_binary(jsonFile, 'continuous', 1);

%retrieve trial onset/offset
E = load_open_ephys_binary(jsonFile, 'events', 1);

startFr = min(C.Timestamps(1), E.Timestamps(1));
endFr = max(C.Timestamps(end), E.Timestamps(end));
fraxis = double(startFr:endFr);
%fraxis = double(C.Timestamps); %this is frame number not time

srate = E.Header.sample_rate;
%fr0 = fraxis(1);
%t0 = fr0/srate;
taxis = 1/srate*fraxis; %time from oephys acquisition start [s]

%does C.Timestamps and E.Timestamps have a common clock?
%if so, E must starts recording earlier than C??

