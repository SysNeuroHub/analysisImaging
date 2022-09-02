function [taxis, trace, ON_frame, OFF_frame] = getDigitalTrace(jsonFile,tgtCh, index)
%[taxis, trace, ON_frame, OFF_frame] = getDigitalTrace(jsonFile,tgtCh, index)

%taxis:
%trace: 
%ON_frame:
%OFF_frame:
%
% TODO
% make this function faster

E = load_open_ephys_binary(jsonFile, 'events', index);

%Look for folder
json=jsondecode(fileread(jsonFile));
header=json.continuous(1);

f=java.io.File(jsonFile);
if (~f.isAbsolute())
    f=java.io.File(fullfile(pwd,jsonFile));
end
f=java.io.File(f.getParentFile(),fullfile('continuous', header.folder_name));
if(~f.exists())
    error('Data folder not found');
end

folder = char(f.getCanonicalPath());

timestamps = readNPY(fullfile(folder,'timestamps.npy'));
taxis = double(timestamps) / double(header.sample_rate);
fr0 = timestamps(1);
nframes = length(taxis);

theseEvents = find(E.ChannelIndex == tgtCh);
tr_ev =E.Data(theseEvents);
tr_fr = E.Timestamps(theseEvents)-fr0+1; %added 1 - is this correct?
ngidx = find(diff(tr_fr)<100);%hack may not be needed
goodidx = setxor(1:length(tr_ev), ngidx);
tr_ev = tr_ev(goodidx);
tr_fr = tr_fr(goodidx);



ON_frame = tr_fr(tr_ev>0); %frame number of trial onset
%ON_times = timestamps(ON_frame); %time of trial offset
OFF_frame = tr_fr(tr_ev<0); %frame number of trial onset
%OFF_times = timestamps(OFF_frame); %time of trial offset

%hack to omit events outside of the timestamps
ON_frame = ON_frame(ON_frame <= nframes);
OFF_frame = OFF_frame(OFF_frame <= nframes);

DItr = zeros(nframes,1);
DItr(ON_frame) = 1;
DItr(OFF_frame) = -1;
trace = cumsum(DItr);
