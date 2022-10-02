function [onTimes, offTimes] = getDigitalTimesOE(jsonFile, digitalCh)
%[onTimes, offTimes] = getDigitalTimesOE(jsonFile, digitalCh)
%retrieve onset and offset of digital signal in digitalCh of jsonFile
%created from getStrobeTimes and acuteDataBasic.getDigitalTimes

if nargin < 2
    digitalCh = 1;
end

%retrieve trial onset/offset
E = load_open_ephys_binary(jsonFile, 'events', 1);
[taxis,fraxis] = getTAxisOE(jsonFile);

try
    theseEvents = find(E.ChannelIndex == digitalCh);
catch err
    onTimes = []; offTimes = [];
    return;
end
%tr_ev =E.Data(theseEvents);
binChar = dec2bin(E.FullWords(theseEvents));
binLogical = binChar=='1';
if ~isempty(binLogical)
    binLogical = fliplr(binLogical);
    
    binLogical_thisCh = binLogical(:,digitalCh);
    
    
    tr_fr = E.Timestamps(theseEvents);% - fr0;
    % ngidx = find(diff(tr_fr)<100);
    % goodidx = setxor(1:length(tr_fr), ngidx);
    % tr_ev = tr_ev(goodidx);
    % tr_fr = tr_fr(goodidx);
    
    
    OnFrames = tr_fr(binLogical_thisCh==1); %frame number of trial onset
    onTimes = taxis(OnFrames-fraxis(1)+1)'; %time of trial onset
    
    OffFrames = tr_fr(binLogical_thisCh==0); %frame number of trial offset ... needs amendment
    offTimes = taxis(OffFrames-fraxis(1)+1)'; %time of trial offset ... needs amendment
else
    onTimes  = [];
    offTimes = [];
end