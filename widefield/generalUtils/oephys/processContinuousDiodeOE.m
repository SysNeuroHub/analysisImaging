function [tOn, tOff, pd3, pdOn, t_ori, t_ds] = processContinuousDiodeOE(jsonFile, channel, sampTo, frScreen, onThresh)
% Process the analog photodiode signal.
% created from marmodata\+marmodata\+cfgs\+acute\neuropix.m
%
% Optional name-value arguments:
%
%   channel - which channel the PD is on (default: 1)
%   sampTo - rate to downsample to (default: 1000 Hz)
%   frScreen - frame rate of the screen (default: 120 Hz)
%   onThresh - threshold for photodiode "on" (default: [], estimated automatically)
%
% Outputs:
%
%   tOn        photodiode on times [s]?
%   tOff:       photodiode off times
%   pd3:        continous, downsampled, filtered trace
%   pdOn:       thresholded, version of pd3.
%   t0:         time stamps before downsampling [s]?
%
% Note: from the original marmodata.neuropix class

%       p = inputParser();
%       p.KeepUnmatched = true;
%       p.addParameter('channel', o.diode, @(x) validateattributes(x,{'double'},{'scalar','nonempty','positive'}));
%       p.addParameter('sampTo', o.sampTo, @(x) validateattributes(x,{'double'},{'scalar','nonempty','positive'})); %hz
%       p.addParameter('frScreen', o.frScreen, @(x) validateattributes(x,{'double'},{'scalar','nonempty','positive'})); %hz
%       p.addParameter('onThresh', o.onThresh, @(x) validateattributes(x,{'double'},{'scalar','nonempty'}));
%
%       p.parse(varargin{:});
%       args = p.Results;
%       ch = args.channel; % photodiode channel is 1
%       frDown = args.sampTo; %hz
%       frScreen = args.frScreen; %fps

if nargin < 5
    onThresh = [];
end

if nargin < 4
    frScreen = 120;
end

if nargin < 3
    sampTo = 1e3;
end
%[pd,t0,frSamp] = o.loadADC(args.channel);
C = load_open_ephys_binary(jsonFile, 'continuous', 1);
pd = C.Data(channel,:);
frSamp = C.Header.sample_rate;
t_ori = double(C.Timestamps)/frSamp; % seconds


guessThresh = onThresh;
if isempty(guessThresh)
    guessThresh = mean(pd);
end

% downsample the continuous signal
%t_ds = decimate(t_ori, floor(frSamp/sampTo));
t_ds = t_ori(1):1/sampTo:t_ori(end);
pd2 = interp1(t_ori,pd, t_ds);

% median filter over 1 frame time
pd3 = movmedian(pd2,floor(1000/frScreen));
pdOn = pd3>guessThresh;

deltaPD = diff(pdOn);
tOnIdx = find(deltaPD==1);
tOffIdx = find(deltaPD==-1);

%[~, tref] = getTAxisOE(jsonFile);
tOn = t_ds(tOnIdx);
tOff = t_ds(tOffIdx);
