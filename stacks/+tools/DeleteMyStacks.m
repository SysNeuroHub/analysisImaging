function DeleteMyStacks(StackDir, FileName)

% DeleteMyStacks deletes generic Stack, saved by SaveMyStacks, from disk
% If the stacks have not yet been saved, it returns error.
%
% ThisStack = DeleteMyStacks([],FileName) loads a stack from the current
% directory
%
% See also: SaveMyStacks LoadMyStacks

% 2015-09-02 DS created from LoadMyStacks

if isempty(StackDir)
    StackDir = pwd;
end



if strcmp(FileName(end-3:end),'.mat')
    FileName = FileName(1:end-4);
end

FileName_mat = sprintf('%s.mat', FileName);
FileName_bin = sprintf('%s.bin', FileName);

try
    % loads something called ThisStack
    if exist(fullfile(StackDir, FileName_bin), 'file') && exist(fullfile(StackDir, FileName_mat), 'file')
        delete(fullfile(StackDir, FileName_bin));
        disp(['Deleted ' fullfile(StackDir, FileName_bin)]);
        
        delete(fullfile(StackDir, FileName_mat));
        disp(['Deleted ' fullfile(StackDir, FileName_mat)]);
    
    elseif exist(fullfile(StackDir, FileName_mat), 'file') % for backward compatibility
        delete(fullfile(StackDir, FileName_mat));
        disp(['Deleted ' fullfile(StackDir, FileName_mat)]);

    else
        disp([FileName_mat ' was not found.']);
    end
    
  
catch err
    disp(['DeleteMyStacks:' err.message]); % added by MC on 2014-05-13
end

