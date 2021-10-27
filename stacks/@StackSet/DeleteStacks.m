function DeleteStacks( S, ServerDir )
% DeleteStacks saves the stacks for an experiment 
%
%   DeleteStacks( ServerDir )
%
% See also: StackSet.LoadStacks, StackSet.SaveStacks

AnimalDir   = fullfile(ServerDir, S.Info.animal);
SeriesDir   = fullfile(AnimalDir, sprintf('%03d', S.Info.iseries));
ExpDir      = fullfile(SeriesDir, sprintf('%03d', S.Info.iexp));

if ischar(S.Info.Stimulus)
    StackDir = fullfile(ExpDir, sprintf('Resize %d %s',...
        round(S.ResizeFactor*100), ...
        S.Info.Stimulus));
else
    StackDir = fullfile(ExpDir, sprintf('Resize %d Stim %03d', ...
        round(S.ResizeFactor*100), ...
        S.Info.Stimulus));
end

if ~exist( StackDir, 'dir')
    warndlg(sprintf('There is no directory %s',StackDir),'Nothing to delete','modal'); 
    exit
end

Choice = questdlg(sprintf('Delete %d files in %s?',S.nConds,StackDir), 'Warning', 'Yes', 'No', 'No');
switch Choice
    case 'No'
        exit;
    otherwise
        for iCond = 1:S.nConds
            fprintf('Loading condition %d of %d from %s\n', iCond, S.nConds, StackDir);
            ThisFileName = sprintf('Stack%04d.mat',iCond);
            delete(fullfile(StackDir, ThisFileName));
        end
        rmdir(StackDir);
end



