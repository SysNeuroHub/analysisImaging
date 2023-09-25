function filtV = hpFilt(V, Fs, cutoffFreq, option)
% filtV = hpFilt(V, Fs, cutoffFreq)
% returns V after temporal filtering (applied to the 2nd dimension) above cutoffFreq
% recommended cutoffFreq = 0.01 (Hz). 
%
% INPUT:
%V: time x space
%
% see also. filtV


if nargin < 4
    option = 0;
end

order = 3;
Wn = cutoffFreq/(Fs/2);
[b,a]=butter(order, Wn, 'high');

% Vtmp = double(V'); %space x time
if option==0
    filtV = single(filtfilt(b,a,double(V')))';
elseif option==1
    filtV = single(zeros(size(V)));
    for ispace = 1:size(V,1)
        filtV(:,ispace) = single(filtfilt(b,a,double(V(ispace,:)')))';
    end
end
    