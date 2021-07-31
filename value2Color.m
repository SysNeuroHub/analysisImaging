function colors = value2Color(values, clim, mapname)
%12/11/20 created from stimID2color

if nargin < 3 
    %colormap(mapname);
    mapname = 'jet';
end

valSize = size(values);

%normSequence = round(255*(value/maxstimID)+1); %stim ID normalized to [1 256]
normSequence = round(255/diff(clim)*(values - clim(1))+1);
normSequence(normSequence>256) = 256;
normSequence(normSequence<1) = 1;
colorfunc = str2func(mapname);
rectColormap = colorfunc(256);
colors = rectColormap(normSequence,:);

colors = reshape(colors, [valSize 3]);
