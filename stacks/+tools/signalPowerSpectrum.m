function [ frequency, power, spectroObj ] = signalPowerSpectrum( signal, sampleFreq, mode )
% return powerspectrum of a time trace

%   WALTHER May 2010


switch nargin
    case {0, 1}
        error('SIGNALPOWERSPECTRUM: input error');
    case 2
        mode = 'mss';
end

if ~strcmp(mode, 'mss') && ~strcmp(mode, 'psd')
    mode = 'mss';
end
    
    
% ******  FTT transformation  ***********
    
% choose periodogram to evaluate the FFT over the full sample length
hspectrum = spectrum.periodogram('hamming');  %('Blackman'); %('hamming');
    
% MSS for the power spectrum and PSD for the power density
if strcmp(mode, 'psd')
    spectroObj = psd(hspectrum, signal,'Fs', sampleFreq, 'NFFT', 'auto');
else
    spectroObj = msspectrum(hspectrum, signal, 'FS', sampleFreq, 'NFFT', 'auto');
end

frequency = spectroObj.Frequencies;
power = spectroObj.Data; 

 %{ 
    % Alternative FFT without windowing
    dt = 1/sampleFreq;
    frequency = ( (0:length(signal)-1)' / ( length(signal) * dt ) )';
   
    signalFT = fft(signal, length(signal))';
    power = ( signalFT .* conj(signalFT) / length(signal) )';
    spectroObj = [];
 %} 
end

