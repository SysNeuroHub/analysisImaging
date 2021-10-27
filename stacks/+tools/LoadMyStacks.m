function ThisStack = LoadMyStacks(StackDir, FileName)

% LoadMyStacks loads generic Stack, saved by SaveMyStacks
% If the stacks have not yet been saved, it returns error.
%
% ThisStack = LoadMyStacks([],FileName) loads a stack from the current
% directory
%
% See also: SaveMyStacks

% 2014-04-02 DS created from LoadStacks
% 2014-05-13 MC made it call tools.loadArr instead of Chris's toolbox
% 2014-11-30 DS allowed FileName to have suffix
% 2014-12-01 DS added error message when FileName does not contain stackset

if isempty(StackDir)
    StackDir = pwd;
end



if strcmp(FileName(end-3:end),'.mat')
    FileName = FileName(1:end-4);
end

FileName_mat = sprintf('%s.mat', FileName);
FileName_bin = sprintf('%s.bin', FileName);

errorCount = 1;
while errorCount <= 5 %try loading 5 times until succeeded
    try
        % loads something called ThisStack
        if exist(fullfile(StackDir, FileName_bin), 'file') && exist(fullfile(StackDir, FileName_mat), 'file')
            disp(['Loading ' fullfile(StackDir, FileName_bin)]);
            [arr, ThisStack] = tools.loadArr(fullfile(StackDir, FileName)); %26/3/14 DS
            ThisStack.Values = arr;
        elseif exist(fullfile(StackDir, FileName_mat), 'file') % for backward compatibility
            disp(['Loading ' fullfile(StackDir, FileName_mat)]);
            load(fullfile(StackDir, FileName_mat)); %load something called ThisStack
        else
            disp([FileName_mat ' was not found.']); break;
        end
        
        if exist('ThisStack','var') break; 
        else
            error(['Check FileName: ' FileName ' does not contain StackSet.']);
        end
    catch err
        disp(['LoadMyStacks:' err.message]); % added by MC on 2014-05-13
        errorCount = errorCount + 1;
        pause(0.1);
        if errorCount == 5
            error('Could not read StackSet');
        end
    end
end

