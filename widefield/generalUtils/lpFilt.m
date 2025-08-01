function filtV = lpFilt(V, Fs, cutoffFreq,     filterEachPixel)
% filtV = ipFilt(V, Fs, cutoffFreq, filterEachPixel)
% returns V after temporal filtering (applied to the 2nd dimension) below cutoffFreq
%
% INPUT:
%V: time x space
%
% see also. filtV

if nargin < 4
    filterEachPixel = 0;
end

order = 3;
Wn = cutoffFreq/(Fs/2);
[b,a]=butter(order, Wn, 'low');

if filterEachPixel == 0
    filtV = single(filtfilt(b,a,double(V')))';
elseif filterEachPixel == 1
    disp('computing lpFilt for each pix...')
    filtV = single(zeros(size(V)));
    for ispace = 1:size(V,1)
        filtV(ispace,:) = single(filtfilt(b,a,double(V(ispace,:)')))';
    end
end