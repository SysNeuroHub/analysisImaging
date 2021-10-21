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

[b100s, a100s] = butter(2, highpassCutoff/(Fs/2), 'high');
[bHeart, aHeart] = butter(2, heartbeatBandStop/(Fs/2), 'stop');

dV = detrend(V', 'linear')';
fVHeart = filter(bHeart,aHeart,dV,[],2);
fV = filter(b100s,a100s,fVHeart,[],2);