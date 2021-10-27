function SaveMyStacks( S, StackDir, LocalTempDir, FileName, doMonitor)
% SaveStacks saves generic stack S to a specified directory as a specified filename
%
%   SaveMyStacks(S, [],[],FileName) saves stack S under the current
%   directory
%
%   SaveMyStacks(S, StackDir,[], FileName) saves stack S under StackDir as
%   FileName.bin (S.Values) and FileName.mat (any other properties)
%
%   SaveMyStacks(S ServerDir, LocalTempDir, FileName ) specifies local temp directory
%
% See also: SaveStacks
%
% 2014-04-02 DS created from SaveStacks.m
% 2014-05-13 MC made it call tools.saveArr instead of Chris's toolbox
% 2015-09-06 DS added 5th input to monitor data size (transient)
% 2016-03-09 DS make directory with random name under temp directory to support parfor

% TO DO: delete temporary files from RAM
% saving of .bin via tools.saveArr can be detected by F-secure antivirus
% use FileRename instead of movefiles, to make it faster?
%
% this function sometimes does not properly work!!!

if nargin < 5
    doMonitor = false;
end

if isempty(StackDir)
    StackDir = pwd;
end

if isempty(LocalTempDir)
    LocalTempDir = 'C:\TEMP';
end

if ~exist(LocalTempDir,'dir')
    LocalTempDir = 'C:\Windows\Temp';
end

%2015/3/4
if strcmp(FileName(end-3:end),'.mat')
    FileName = FileName(1:end-4);
end

if ~exist(StackDir, 'dir') %13/8/15
    mkdir(StackDir);
end

if doMonitor
    oriStackInfo = whos('S');
    oriStackBytes = oriStackInfo.bytes;
    disp(['SaveMyStacks:original bytes:' num2str(oriStackBytes)]);
end

arr = S.Values;
S.Values = [];

tempname = randfname(LocalTempDir, 50);
mkdir(fullfile(LocalTempDir,tempname));
tools.saveArr( fullfile(LocalTempDir,tempname,tempname), arr, S);%this is sometimes detected by F-secure antivirus

if doMonitor
    tempfileInfo = dir([fullfile(LocalTempDir,tempname,tempname),'.bin']);
    tempStackBytes = tempfileInfo.bytes;
    disp(['SaveMyStacks:temp file bytes:' num2str(tempStackBytes)]);
end

errorCount = 1;
while errorCount <= 5 %try transferring 5 times until succeeded. %16/5/14 DS
    try
        movefile( [fullfile(LocalTempDir,tempname,tempname) '.bin'], [fullfile(StackDir,FileName) '.bin']);
        movefile( [fullfile(LocalTempDir,tempname,tempname) '.mat'], [fullfile(StackDir,FileName) '.mat']);
        break;
    catch err
        disp(err.message); % added by MC on 2014-05-13
        errorCount = errorCount + 1;
        pause(errorCount);%1/9/15
        if errorCount > 6
            error('Could not move StackSet');
        end
    end
end

if doMonitor
    try
        finalfileInfo = dir([fullfile(StackDir,FileName) '.bin']);
        finalStackBytes = finalfileInfo.bytes;
        disp(['SaveMyStacks:final bytes:' num2str(finalStackBytes)]);
        
        if (oriStackBytes/finalStackBytes > 1.5) || (oriStackBytes/finalStackBytes < 0.5)
            disp('SaveMyStacks: re-saving stacks...');
            SaveMyStacks( S, StackDir, LocalTempDir, FileName);
        end
    catch err
        disp('SaveMyStacks: re-saving stacks...');
        SaveMyStacks( S, StackDir, LocalTempDir, FileName);
    end
end

rmdir(fullfile(LocalTempDir,tempname),'s'); %delete temp directory

end

function fname = randfname(dirname, sLength)
%generate unique random sequence for filename in tempdir, so
%SaveStacks can be called at the same time in a single host

s = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789';

%find number of random characters to choose from
numRands = length(s)-1;

rng('shuffle','twister');%DS 2015/3/29

%generate random string
fname = [s( round(rand(1,sLength)*numRands) + 1)];
while exist(fullfile(dirname, fname), 'file')
    fname = [s( round(rand(1,sLength)*numRands) +1)];
end
end