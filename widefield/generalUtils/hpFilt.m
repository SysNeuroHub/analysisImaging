function filtV = hpFilt(V, Fs, cutoffFreq)
% filtV = hpFilt(V, Fs, cutoffFreq)
% returns V after temporal filtering (applied to the 2nd dimension) above cutoffFreq
% recommended cutoffFreq = 0.01 (Hz). 

order = 3;
Wn = cutoffFreq/(Fs/2);
[b,a]=butter(order, Wn, 'high');

filtV = single(filtfilt(b,a,double(V')))';