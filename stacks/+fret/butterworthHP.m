function [Hd, cutoffFreq] = butterworthHP(Fsampling, Fstop, Fpass, Astop)
%Hd = butterworthHP(Fsampling, Fstop, Fpass, Astop) returns a discrete-time filter object.
% Fstop = 0.3;  % Stopband Frequency
% Fpass = 0.5;  % Passband Frequency
% Astop = 40;   % Stopband Attenuation (dB)
%
% M-File generated by MATLAB(R) 7.9 and the Signal Processing Toolbox 6.12.
% 2013.12.17 Generated by DS
% 2018.1.4 added 2nd output cutoffFreq

% Butterworth Highpass filter designed using the BUTTER function.

% All frequency values are in Hz.
Fs = Fsampling;  % Sampling Frequency

Apass = 1;    % Passband Ripple (dB)

% Calculate the order from the parameters using BUTTORD.
[N,Fc] = buttord(Fpass/(Fs/2), Fstop/(Fs/2), Apass, Astop);

% Calculate the zpk values using the BUTTER function.
[z,p,k] = butter(N, Fc, 'high');
cutoffFreq = Fc*Fs/2;

% To avoid round-off errors, do not use the transfer function.  Instead
% get the zpk representation and convert it to second-order sections.
[sos_var,g] = zp2sos(z, p, k);
%freqz(sos_var,2^10,Fs);%frequency response
Hd          = dfilt.df2sos(sos_var, g);

% [EOF]