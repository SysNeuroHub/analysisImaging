function filtV = lpFilt(V, Fs, cutoffFreq)
% filtV = hpFilt(V, Fs, cutoffFreq)
% returns V after temporal filtering above cutoffFreq
% recommended cutoffFreq = 0.01 (Hz). 

order = 3;
Wn = cutoffFreq/(Fs/2);
[b,a]=butter(order, Wn, 'low');

try
    filtV = single(filtfilt(b,a,double(V')))';
catch err
    disp('computing lpFilt for each pix...')
    filtV = single(zeros(size(V)));
    for ispace = 1:size(V,2)
        filtV(:,ispace) = single(filtfilt(b,a,double(V(:,ispace)')))';
    end
end