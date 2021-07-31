function stimIDColor = stimID2color(stimID, maxstimID, mapname)
%stimIDColor = stimID2color(stimID, maxstimID, mapname) [ R G B ]
%15/6/20 created 

if nargin < 3 
    %colormap(mapname);
    mapname = 'jet';
end

normSequence = round(255*(stimID/maxstimID)+1); %stim ID normalized to [1 256]
colorfunc = str2func(mapname);
%rectColormap = colormap;
rectColormap = colorfunc(256);
stimIDColor = rectColormap(normSequence,:);
