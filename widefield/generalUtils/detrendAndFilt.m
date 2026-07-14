function fV = detrendAndFilt(V, Fs, highpassCutoff, heartbeatBandStop)
% fV = detrendAndFilt(V, Fs, highpassCutoff, heartbeatBandStop)
% removes baseline.
% output will be dV rather than V

if nargin < 3
    highpassCutoff = 0.2;%4/12/17 DS 0.01; % Hz
end
if nargin < 4
    heartbeatBandStop = [9 14];
end

dV = detrend(V', 'linear')'; %this will make the mean over time 0

if sum(~isinf(heartbeatBandStop))==2
    [bHeart, aHeart] = butter(2, heartbeatBandStop/(Fs/2), 'stop');
    fVHeart = filter(bHeart,aHeart,dV,[],2);
else
    fVHeart = dV;
end
if ~isinf(highpassCutoff)
    [b100s, a100s] = butter(2, highpassCutoff/(Fs/2), 'high');
    fV = filter(b100s,a100s,fVHeart,[],2);
else
    fV = fVHeart;
end