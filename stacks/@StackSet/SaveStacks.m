function [StackDir, ThisFileName] = SaveStacks( S, ServerDir, LocalTempDir, FileString, suffix)
% StackSet.SaveStacks(ServerDir) 
% saves StackSet to file using tools.saveArr
% Data is stored at:
% //ServerDir/animal/iseries/iexp/
% 
% StackSet.SaveStacks(ServerDir, LocalTempDir, FileString) 
%  allows to specify a local temp directory (DEFAULT:  'C:\TEMP').
% 
% StackSet.SaveStacks(ServerDir, LocalTemp, FileString) 
% changes directory to save as:
% //ServerDir/[animal FileString]/iseries/iexp/
% 
% Tips: StackSet data is separately stored into two files: OOO.bin (StackSet.Values) and OOO.mat (other StackSet properties).
% 
% See also: tools.saveArr, StackSet.SaveStacks, StackSet.GetOneCondition%
% 2014.3.27 DS S.Values is in single precision
% 2015.1.27 DS added 5th input (suffix)
% 2015.9.2 DS added 1st and 2nd outputs

S.Values = single(S.Values); %27/3/14 DS

if nargin < 5 %27/1/15
    suffix = '';
end

if nargin < 4 %13.10.18. DS. maybe not good idea to have the 4th input
    FileString = '';
end

if nargin < 3
    LocalTempDir = 'C:\TEMP';
end

if isempty(LocalTempDir)
    LocalTempDir = 'C:\TEMP';
end

if ~exist(LocalTempDir,'dir')
    LocalTempDir = 'C:\Windows\Temp';
end

AnimalDir   = fullfile(ServerDir, [S.Info.animal FileString]);
if isstr(S.Info.iseries)
SeriesDir   = fullfile(AnimalDir, sprintf('%s', S.Info.iseries));
else
SeriesDir   = fullfile(AnimalDir, sprintf('%03d', S.Info.iseries));
end
ExpDir      = fullfile(SeriesDir, sprintf('%03d', S.Info.iexp));

% if ischar(S.Info.Stimulus)
%     StackDir = fullfile(ExpDir, sprintf('Resize %d %s',...
%         round(S.ResizeFactor*100), ...
%         S.Info.Stimulus));

if strcmp(S.Info.Stimulus, 'All') && strcmp(S.Info.RepeatList, 'Average') %all stim, avg repeats
    LoadOption = 1;
elseif ~strcmp(S.Info.Stimulus, 'All') && strcmp(S.Info.RepeatList, 'Average')%single stim, avg repeats
    LoadOption = 2;
elseif ~strcmp(S.Info.Stimulus, 'All') && ~isempty(S.Info.RepeatList) %single repeat, single stim
    LoadOption = 3;
elseif strcmp(S.Info.Stimulus, 'All') &&  ~isempty(S.Info.RepeatList) %all stim, single repeat
    LoadOption = 4;
end

switch LoadOption
    case 1
        StackDir = fullfile(ExpDir, sprintf('Resize %d Stim All Repeat Average',...
            round(S.ResizeFactor*100)));
    case 2
        StackDir = fullfile(ExpDir, sprintf('Resize %d Stim %03d Repeat Average',...
            round(S.ResizeFactor*100), ...
            S.Info.Stimulus));
    case 3
        StackDir = fullfile(ExpDir, sprintf('Resize %d Stim %03d Repeat All', ...
            round(S.ResizeFactor*100), ...
            S.Info.Stimulus));
    case 4
        StackDir = fullfile(ExpDir, sprintf('Resize %d Stim All Repeat %03d', ...
            round(S.ResizeFactor*100),S.Info.RepeatList)); %21.01.14 AP
end

if ~isempty(suffix)
    StackDir = [StackDir ' ' suffix];
end

if ~exist(AnimalDir, 'dir'), mkdir(AnimalDir); end
if ~exist(SeriesDir, 'dir'), mkdir(SeriesDir); end
if ~exist(   ExpDir, 'dir'), mkdir(   ExpDir); end
if ~exist( StackDir, 'dir'), mkdir( StackDir); end

switch LoadOption
    case {1, 2, 4}
        %option1,2: suffix from S.Info.Stimulus?
        %option4: suffix from S.Info.RepeatList?
        fprintf('Saving %d files in %s\n',S.nConds,StackDir);
        for iCond = 1:S.nConds
            ThisStack = S.GetOneCondition(iCond); % ignore the warning
            fprintf('Saving condition %d of %d\n',iCond,S.nConds);
            ThisFileName = sprintf('Stack%04d',iCond);
        end
    case 3
        fprintf('Saving a file in %s\n',StackDir);
        ThisStack = S;
        %ThisStack = S.GetOneCondition(iCond); % ignore the warning
        if ~isstr(ThisStack.Info.RepeatList)
            fprintf('Saving condition %d \n', ThisStack.Info.RepeatList);
            %ThisFileName = sprintf('Stack%04d.mat', ThisStack.Info.RepeatList);
            ThisFileName = sprintf('Stack%04d', ThisStack.Info.RepeatList);
        else %for saving median across trials. %13.12.23 - DS
            fprintf('Saving condition %s \n', ThisStack.Info.RepeatList);
            %ThisFileName = sprintf('Stack%s.mat', ThisStack.Info.RepeatList);
            ThisFileName = sprintf('Stack%s', ThisStack.Info.RepeatList);
        end
end

%             save( fullfile(LocalTempDir,'TempStack.mat'), 'ThisStack','-v7.3');
%             movefile( fullfile(LocalTempDir,'TempStack.mat'), fullfile(StackDir,ThisFileName) );
arr = ThisStack.Values;
ThisStack.Values = [];

tempname = randfname(LocalTempDir, 50);
tools.saveArr( fullfile(LocalTempDir,tempname), arr, ThisStack);
%add catch-error statement
errorCount = 1;
while errorCount <= 5 %try saving 5 times until succeeded. DS added on 24/4/14
    try
        movefile( [fullfile(LocalTempDir,tempname) '.bin'], [fullfile(StackDir,ThisFileName) '.bin']);
        movefile( [fullfile(LocalTempDir,tempname) '.mat'], [fullfile(StackDir,ThisFileName) '.mat']);
        disp(' ');
        
        break;
        
    catch err
        errorCount = errorCount + 1;
        pause(errorCount);%1/9/15
        LocalTempDir
        tempname
        StackDir
        ThisFileName
    end
end
end
%clear ThisStack %2014-1-7 DS


function fname = randfname(dirname, sLength)
%generate unique random sequence for filename in tempdir, so
%SaveStacks can be called at the same time in a single host

s = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789';

%find number of random characters to choose from
numRands = length(s)-1;

rng('shuffle','twister');%DS 2015/3/27

%generate random string
fname = [s( round(rand(1,sLength)*numRands) + 1)];
while exist(fullfile(dirname, fname), 'file')
    fname = [s( round(rand(1,sLength)*numRands) +1)];
end
end