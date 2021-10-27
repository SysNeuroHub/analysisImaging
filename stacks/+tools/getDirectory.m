function [StackDir, ThisFileName] = getDirectory( ServerDir, p, ResizeFactor, iStim, RepeatList, Cam, suffix)
% [StackDir] = getDirectory( ServerDir, p, ResizeFactor, iStim, RepeatList, Cam) 
% returns directory of stackset, stored under ServerDir
% [StackDir, ThisFileName] = getDirectory( ServerDir, p, ResizeFactor, iStim, RepeatList, Cam) 
% also returns name filename with .mat suffix
% [StackDir] = getDirectory( ServerDir, p, ResizeFactor, iStim, RepeatList, Cam, suffix) 
% returns directory name, with specified suffix.
% eg. \\zserver3\Data\Stacks\Kahvesi46_ratio\20141019\003\Resize 50 Stim 001 Repeat All nohp 
%
% See also: tools.getServerDirectory (returns directory of raw data)

% only Cam.FileString is used...
% 2013-12-31 DS created from LoadStacks
% 2014-12-07 DS added "suffix" input
global DIRS

if nargin < 7
    suffix = '';
end
 
if nargin < 6
    Cam.FileString=[];
end

if iscell(Cam.FileString) %2015/3/27
    Cam.FileString = cell2mat(Cam.FileString);
end

if nargin < 5
    RepeatList = [];
end

if nargin < 4
    iStim = [];
    RepeatList = [];
end

if nargin < 3
    ResizeFactor = 1;
end


if isempty(iStim) && isempty(RepeatList)
    LoadOption = 1; %all stimuli, avg repeats
elseif ~isempty(iStim) && isempty(RepeatList)
    LoadOption = 2; % single stimulus, average repeats
elseif ~isempty(iStim) && ~isempty(RepeatList)
    LoadOption = 3; % single stimulus, specified repeat(s?)
else
    error('cannot recognize load option')
end

AnimalDir   = fullfile(ServerDir, [p.animal Cam.FileString]); % 13.10.18 DS
if isnumeric(p.iseries)
    SeriesDir   = fullfile(AnimalDir, sprintf('%03d', p.iseries));
else
    SeriesDir   = fullfile(AnimalDir, sprintf('%s',p.iseries));
%     SeriesDir   = fullfile(AnimalDir, sprintf('%s', ['20' num2str(p.iseries(3:4)) '-' num2str(p.iseries(5:6)) ...
%         '-' num2str(p.iseries(7:8))]));
end
ExpDir      = fullfile(SeriesDir, sprintf('%03d', p.iexp));

switch LoadOption
    case 1
        StackDir = fullfile(ExpDir, sprintf('Resize %d Stim All Repeat Average',...
            round(ResizeFactor*100)));
        ThisFileName = sprintf('Stack%04d.mat',iStim);
    case 2
        StackDir = fullfile(ExpDir, sprintf('Resize %d Stim %03d Repeat Average', ...
            round(ResizeFactor*100), iStim) );
        ThisFileName = sprintf('Stack%04d.mat',1);
    case 3
        StackDir = fullfile(ExpDir, sprintf('Resize %d Stim %03d Repeat All', ...
            round(ResizeFactor*100), iStim) );
        ThisFileName = sprintf('Stack%04d.mat',RepeatList); %is this correct?
end

if ~isempty(suffix)
    StackDir = [StackDir ' ' suffix];
end