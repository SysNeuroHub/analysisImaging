function filtV = lpFilt(V, Fs, cutoffFreq)
% filtV = hpFilt(V, Fs, cutoffFreq)
% returns V after temporal filtering above cutoffFreq
% recommended cutoffFreq = 0.01 (Hz). 

order = 3;
Wn = cutoffFreq/(Fs/2);
[b,a]=butter(order, Wn, 'low');

filtV = single(filtfilt(b,a,double(V')))';