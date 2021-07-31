function [ times, lfp ] = grabLFP( expt, expr, frames, TSDir)
% [ times ] = grabFrameTimes( expt, expr, frames, TSDir)
% returns times of TTL signals going from LOW to HIGH
% for 2p frame time (default), return times of each plane as a cell (planes x 1)
% for other TSFieldNames ,return an array
%23/3/20 created from grabFrameTimes

subject = expt.subject;
expDate = expt.expDate;

if nargin < 4
    TSDir = 'D:\thorimagedata\'; %where ThorSync data is saved or copied
end

if nargin < 3
    frames = 'all';
end

%timelineFile = strcat(expDate,'_',num2str(expr),'_',subject,'_','Timeline.mat');
fnam = fullfile(TSDir, subject, expDate, num2str(expr), 'Episode001.h5');
[syncDataOut] = LoadSyncEpisodeFunction(fnam);


times = syncDataOut.time;


if isfield(syncDataOut, 'npi')
   lfp = syncDataOut.npi;
else
    warning('LFP is not recorded');
end


