function V = filtV(V, Fs, cutoffFreq, lpFreq, filterEachPix)
% V = filtV(V, Fs, cutoffFreq, lpFreq)
% apply low cutt and high pass filters on 2nd dimension of V
%
% see also. hpFilt, lpFilt

if nargin < 3 
    lpFreq = [];
end
if nargin < 4
    filterEachPix = 0;
end

if ~isempty(cutoffFreq)
    meanV = mean(V,2);
    V =  hpFilt(V-meanV, Fs, cutoffFreq, filterEachPix); %cutoffFreq
    disp(['high pass filetered at ' num2str(cutoffFreq) '[Hz]'])
    V = V + meanV;
end
if ~isempty(lpFreq)
    meanV = mean(V,2);%must be temporal avg
    V =  lpFilt(V-meanV, Fs, lpFreq, filterEachPix); %low-pass Freq
    disp(['low pass filetered at ' num2str(lpFreq) '[Hz]'])
    V = V + meanV;
end
