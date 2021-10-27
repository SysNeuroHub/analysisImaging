function [ slaveCorrection, masterCorrection ] = slaveMasterAdjustFactor2D...
    ( slave, master, Fpass, Fsampling, filterSize, firstFrame, lastFrame)
% returns noise-gain equalization factors for each pixel from slave and master sequences
% [ slaveCorrection, masterCorrection ] = slaveMasterAdjustFactor2D( slave,
% master, Fpass, Fsampling, filterSize)
% returns slave and master factors 
%
% Inputs:
% slave, master: x-y-t tensors of the raw signal (F instead of dF/F). slave
% and master needs to be the same size
% Fpass: range of frequency to calculate the gain-equalization [min max] [Hz]
% Fsampling: sampling frequency of the slave and master tensors [Hz]
% filtersize: spatial filter size in pixel
%
% [...] = slaveMasterAdjustFactor2D( slave, master, Fpass, Fsampling,
% filterSize, firstFrame, lastFrame)
% lets you specify period of slave/master tensors to calculate the
% equalization factors.
%
% See also: fret.makeratio


% WALTHER Mar 2011
%
% TO DO: remove edge of master image so that the spatial filtering does not
% give an artifact

if nargin <= 3 || ~isequal(size(Fpass), [1,2]),...
        error('SEQUENCE ADJUST: wrong arguments');
end

if ~isequal(size(slave),size(master))
    error('master and slave sizes do not match');
end

switch nargin
    case {0, 1,2,3}
        error('SEQUENCE ADJUST: input missing');
    case 4
         firstFrame = 1;
        lastFrame = inf;
        filterSize = 1;
    case 5
         firstFrame = 1;
        lastFrame = inf;
    case 6
        lastFrame = inf;
end

frames = size(master,3);

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


%hack to eliminate boundary effect 2/5/14
% masterTh = prctile(master(:),1);
% master(master < masterTh) = masterTh;

%these offset values should be exactly the same to sequenceBaseImage in
%seququenceAdjust2D
masterOffset = squeeze(mean(master(:,:,firstFrame:lastFrame), 3));%modified on 6.21.2012
slaveOffset = squeeze(mean(slave(:,:,firstFrame:lastFrame), 3));%modified on 6.21.2012

[dum,MASTEROFFSET] = meshgrid(ones(size(master,3),1),masterOffset(:));
MASTEROFFSET = reshape(MASTEROFFSET,size(masterOffset,1), size(masterOffset,2), size(master,3));
master = master - MASTEROFFSET;
clear MASTEROFFSET;

[dum,SLAVEOFFSET] = meshgrid(ones(size(slave,3),1),slaveOffset(:));
SLAVEOFFSET = reshape(SLAVEOFFSET,size(slaveOffset,1), size(slaveOffset,2), size(slave,3));
slave = slave - SLAVEOFFSET;
clear SLAVEOFFSET;

%spatial filtering is absolutely necessary! 
%very slow median filter did not really prevent overshoot around bessels.
if size(master,1) ~= 1 || size(master,2) ~= 1
    disp('slaveMasterAdjustFactor2D: spatial filtering..')
    master = tools.averageFilter2(master, filterSize, 'average');%should replace with more generic filter
    slave = tools.averageFilter2(slave, filterSize, 'average');
end

%temporal filtering.
master = filter( Hd, master, 3 );%initial frames can be NAN...why?
slave = filter( Hd, slave, 3 );


masterCorrection = single( 0.5 * ( 1 + ( masterOffset .* nanstd(slave(:,:,firstFrame:lastFrame),0,3) ) ...
    ./ ( slaveOffset .* nanstd(master(:,:,firstFrame:lastFrame),0,3) ) ) );

slaveCorrection = single( 0.5 * ( 1 + ( slaveOffset .* nanstd(master(:,:,firstFrame:lastFrame),0,3) ) ...
    ./ ( masterOffset .* nanstd(slave(:,:,firstFrame:lastFrame),0,3) ) ) );


end


