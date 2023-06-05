function V = filtV(V, Fs, cutoffFreq, lpFreq)
% V = filtV(V, Fs, cutoffFreq, lpFreq)

if ~isempty(cutoffFreq)
    meanV = mean(V,2);
    V =  hpFilt(V-meanV, Fs, cutoffFreq); %cutoffFreq
    disp(['high pass filetered at ' num2str(cutoffFreq) '[Hz]'])
    V = V + meanV;
end
if ~isempty(lpFreq)
    meanV = mean(V,2);%must be temporal avg
    V =  lpFilt(V-meanV, Fs, lpFreq); %low-pass Freq
    disp(['low pass filetered at ' num2str(lpFreq) '[Hz]'])
    V = V + meanV;
end
