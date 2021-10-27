function [ slaveCorrection, masterCorrection ] = slaveMasterAdjustFactor( slave, master, Fpass, Fsampling, mask, firstBaseFrame, lastBaseFrame, firstFrame, lastFrame)
%UNTITLED Summary of this function goes here

% WALTHER Mar 2011
%   Detailed explanation goes here

if nargin <= 4 || ~isequal(size(Fpass), [1,2]) || ~isequal(size(slave),size(master)) || ~isequal(size(master(:,:,1)),size(mask)),...
       error('SEQUENCE ADJUST: wrong arguments');
end

switch nargin
    case {0, 1,2,3,4}
        error('SEQUENCE ADJUST: input missing');
    case 5
        firstBaseFrame = 1;
        lastBaseFrame = inf;
        firstFrame = 1;
        lastFrame = inf;
    case 6
        lastBaseFrame = inf;
        firstFrame = 1;
        lastFrame = inf;
    case 7
        firstFrame = 1;
        lastFrame = inf;
    case 8
        lastFrame = inf;
end

frames = size(master,3);

baseFrame = [firstBaseFrame; lastBaseFrame];
baseFrame( baseFrame > frames ) = frames;
baseFrame( baseFrame < 1 ) = 1;
baseFrame = sort( baseFrame );
firstBaseFrame = baseFrame(1);
lastBaseFrame = baseFrame(2);

limitFrame = [firstFrame; lastFrame];
limitFrame( limitFrame > frames ) = frames;
limitFrame( limitFrame < 1 ) = 1;
limitFrame = sort( limitFrame );
firstFrame = limitFrame(1);
lastFrame = limitFrame(2);


% Design bandpass filter
Fs = Fsampling;  % Sampling Frequency

Fstop1 = Fpass(1)-1;               % First Stopband Frequency
Fpass1 = Fpass(1);               % First Passband Frequency
Fpass2 = Fpass(2);              % Second Passband Frequency
Fstop2 = Fpass(2)+1;              % Second Stopband Frequency
Dstop1 = 0.01;             % First Stopband Attenuation
Dpass  = 0.1;               % Passband Ripple
Dstop2 = 0.01;             % Second Stopband Attenuation
dens   = 20;              % Density Factor

% Calculate the order from the parameters using FIRPMORD.
[N, Fo, Ao, W] = firpmord([Fstop1 Fpass1 Fpass2 Fstop2]/(Fs/2), [0 1 ...
                          0], [Dstop1 Dpass Dstop2]);

% Calculate the coefficients using the FIRPM function.
b  = firpm(N, Fo, Ao, W, {dens});
Hd = dfilt.dffir(b);


% adjust slave sequence

signal1 = zeros(frames,1);
signal2 = zeros(frames,1);
                
for k=1:frames
    dummy = double(master(:,:,k));
    signal1(k) = mean(dummy(mask));
    dummy = double(slave(:,:,k));
    signal2(k) = mean(dummy(mask));
end

signal1Offset = signal1(firstBaseFrame);
signal2Offset = signal2(firstBaseFrame);

%{
signal1Offset = mean(signal1(firstBaseFrame:firstBaseFrame));
signal2Offset = mean(signal2(firstBaseFrame:lastBaseFrame));
%}

signal1 = signal1(firstBaseFrame:lastBaseFrame) - repmat(signal1Offset,lastBaseFrame-firstBaseFrame+1,1);
signal2 = signal2(firstBaseFrame:lastBaseFrame) - repmat(signal2Offset,lastBaseFrame-firstBaseFrame+1,1); 

signal1 = filter( Hd, signal1 );
signal2 = filter( Hd, signal2 );

masterCorrection = single( 0.5 * ( 1 + ( signal1Offset * std(signal2) ) / ( signal2Offset * std(signal1) ) ) );

slaveCorrection = single( 0.5 * ( 1 + ( signal2Offset * std(signal1) ) / ( signal1Offset * std(signal2) ) ) );

end


