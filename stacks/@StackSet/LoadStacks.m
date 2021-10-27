function [S, registerInfo, ratioInfo] = LoadStacks( ServerDir, p, ResizeFactor, ...
    iStim, RepeatList, Cam, varargin )
% 
% LoadStacks loads the stacks for an experiment.
% 
% If the stacks have already been saved, it loads those. Otherwise, it
% loads them from the raw data.
% 
%   S = LoadStacks( ServerDir, protocol ) loads all stimuli for the
%   experiment indicated in protocol (as obtained by LoadAllStimuliAverageRepeats). 
%   Data are averaged across repeats.
% 
%   S = LoadStacks( ServerDir, protocol, ResizeFactor ) lets you specify
%   the ResizeFactor (DEFAULT: 1)
% 
%   S = LoadStacks( ServerDir, protocol, ResizeFactor, iStim ) loads the
%   averaged repeat for stimulus iStim (DEFAULT: [])
% 
%   S = LoadStacks( ServerDir, protocol, ResizeFactor, iStim, RepeatList ) loads the
%   specified individual repeat for stimulus iStim (DEFAULT: []). nFrames
%   can be different across repeats (currently RepeatList must be a scalar. see LoadOneStimulus.m)
% 
%   S = LoadStacks( ServerDir, protocol, ResizeFactor, [], RepeatList ) loads the
%   specified individual repeats for all stimuli
% 
%   S = LoadStacks( ServerDir, protocol, ResizeFactor, iStim, RepeatList, Cam )
%   specifies the imager object. Useful when creating the stacks from raw
%   data. If Cam.FileString = '_ratio', load ratiometric signal
% 
%   S = LoadStacks( __, Name, Value )
%   specifies additional information to make stackset using one or more
%   Name, Value pair arguments. Name, Value pair settings apply to all the
%   stacks loaded.
%   Name-Value pair arguments
%   'ratioInfo' -- information needed to compute ratiometric
%   ‘registerInfo’ -- information needed for image registration (currently implemented only for ratiometric signal)
%   'newSave' -- does not save stackset if false (default:true). added on
%   25.5.14 DS
%   'suffix' -- suffix for directory name. added on 27.1.15 DS
% 
% Tips: LoadStacks loads a StackSet for an experiment from file using loadArr. 
% OOO.bin and OOO.mat files are integrated as “ThisStack”. You can read other 
% properties (eg. timeVec) of stackset other than .Values by just loading the .mat very quickly.
% 
% See also: tools.getDirectory, tools.Imager, StackSet.SaveStacks, 
% StackSet.LoadOneStimulus, StackSet.LoadOneStimulusAverageRepeats, 
% StackSet.LoadAllStimuliAverageRepeats, StackSet.LoadAllStimuliSelectedRepeats

%TO DO: in option 3, enable RepeatList to be a vector

% 2014-12-08 DS flipped 2nd and 3rd output
% 2015-01-27 DS added "suffix" as a varargin input (currently only for ratio)
% 2015-08-14 DS in case of ratiometric, save ratio and sum stacks in separate
% directories

global DIRS


%registerImage = [];
ratioInfo = [];
newSave = true;
if ~isempty(varargin)
    
    vidx = [];
    for vv = 1 : length(varargin)
        if isstr(varargin{vv})
            vidx = [vidx vv];
        end
    end
    
    %information for image registration in varargin
    for vv = vidx
        if any(strfind(varargin{vv}, 'registerInfo'))
            registerInfo = varargin{vv+1};
            break
        end
    end
    
    %information for ratiometric in varargin
    for vv = vidx
        if any(strfind(varargin{vv}, 'ratioInfo'))
            ratioInfo = varargin{vv+1};
            break
        end
    end
    
    %whether to save newly crasted stackset
    for vv = vidx
        if any(strfind(varargin{vv}, 'newSave'))
            newSave = varargin{vv+1};
            break
        end
    end
    
    % any suffix for directory name. 27/1/15 
    for vv = vidx
        if any(strfind(varargin{vv}, 'suffix'))
            suffix = varargin{vv+1};
            break
        end
    end
end

if ~exist('suffix','var')
    suffix = '';
end

if nargin < 6
    FileString = [];
    Cam = [];
else
    FileString = Cam.FileString; %14.01.14 AP
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
elseif isempty(iStim) && ~isempty(RepeatList) % 14.01.21 AP
    LoadOption = 4; % all stimuli, specified repeat
else
    error('cannot recognize load option')
end

[StackDir, ThisFileName] = tools.getDirectory( ServerDir, p, ResizeFactor, iStim, RepeatList, Cam, suffix);

%if LoadOption = 2or3, judge by stackset file.
%if LoadOption = 1or4, judge by directory
if (((LoadOption==2)||(LoadOption==3))&&exist( fullfile(StackDir, ThisFileName), 'file')) || ...
        (((LoadOption==1)||(LoadOption==4))&&exist( StackDir, 'dir')) 
    %% The stacks have already been saved. Load them one by one.    
    fprintf('Found saved stacks for %s series %d experiment %d resize %d in\n', ...
        p.animal,p.iseries,p.iexp, ResizeFactor*100);
    fprintf('%s\n',StackDir);
    
    switch LoadOption
        case 1
            nConds = p.nstim;
            idx = 1:nConds;
            fprintf('Loading %d stimuli (averaged across repeats)',nConds);
        case 2
            nConds = 1;
            idx = 1;
            fprintf('Loading iStim %d (averaged across repeats)', iStim);
        case 3
            nConds = length(RepeatList);
            idx = RepeatList;
            fprintf('Loading %d repeats (for iStim %d)',idx, iStim);%load all repeats??
        case 4 %21.01.14 AP
            nConds = p.nstim;
            idx = 1:nConds;
            fprintf('Loading repeat: %d (for all stimuli)',RepeatList);
    end
    
    
    ii = 1;
    for iCond = idx
        fprintf('.');
        ThisFileName = sprintf('Stack%04d.mat',iCond);
        ThisFileName_bin = sprintf('Stack%04d.bin',iCond);
        
        if exist('ThisStack','var')
            clear ThisStack
        end
        
        errorCount = 1;
        while errorCount <= 5 %try loading 5 times until succeeded
            try
                % loads something called ThisStack
                if exist(fullfile(StackDir, ThisFileName_bin), 'file') && exist(fullfile(StackDir, ThisFileName), 'file')
                   [arr, ThisStack] = tools.loadArr(fullfile(StackDir, sprintf('Stack%04d',iCond))); %26/3/14 DS
                   ThisStack.Values = arr;
                   break;%DS 4/4/14
                elseif exist(fullfile(StackDir, ThisFileName), 'file') % for backward compatibility
                    load(fullfile(StackDir, ThisFileName));
                    break; %DS 4/4/14
                end
                
                if exist('ThisStack','var') break; end
            catch err
                errorCount = errorCount + 1;
                warning('Failed to load StackSet');
                pause(0.1);
            end
        end
        
        if iCond == 1 || nConds == 1;
            S = ThisStack;
            S.nConds = nConds;
            S.Values = nan(S.nRows, S.nCols, S.nFrames, S.nConds, 'single');%27/3/14 DS
        end
        S.Values(:,:,:,ii) = ThisStack.Values(:,:,:,1);%29/11/14
        ii = ii+1;
    end
    fprintf('done.\n');
    
    if nargout == 2 %2014/12/8 modified
        try
            load(fullfile(StackDir, ['info' ThisFileName]));
            registerInfo = info.registerInfo;
        catch err
            warning('Could not read peripheral information');
        end
    end
    if nargout == 3 %2014/12/8 added
        try
            load(fullfile(StackDir, ['info' ThisFileName]));
            ratioInfo = info.ratioInfo;
        catch err
            warning('Could not read peripheral information');
        end
    end
    
    %     fprintf('.');
    %     load(fullfile(StackDir, ThisFileName)); % loads something called ThisStack
    %     S = ThisStack;
    %     S.nConds = nConds;
    %     S.Values = nan(S.nRows, S.nCols, S.nFrames, S.nConds);
    %     %         S.Values(:,:,:,iCond) = ThisStack.Values;
    %     S.Values = ThisStack.Values;
    %     fprintf('done.\n');
    
else
    %% The stacks have not been saved. Load them from raw data.
    %     TheAnswer = questdlg('No stack found. Make a new one from raw data, and save it?', ...
    %         'No stack found', ...
    %         'Yes', 'No', 'Yes');
    display('StackSet was not found. Build StackSet from scratch...')
    TheAnswer = 'Yes';
    switch TheAnswer
        case 'Yes'
            if isempty(Cam)
                Cam = Imager();
            end
            
            switch LoadOption
                case 1
                    S = StackSet.LoadAllStimuliAverageRepeats(Cam, p, ResizeFactor);
                case 2
                    if strcmp(Cam.FileString , '_ratio')
                        S = StackSet.LoadOneStimulusAverageRepeats_ratio(ServerDir, Cam, p, ResizeFactor, iStim, ratioInfo);
                    else
                        S = StackSet.LoadOneStimulusAverageRepeats(Cam, p, ResizeFactor, iStim);
                    end
                case 3
                    if strcmp(Cam.FileString , '_ratio')
                            [S, ratioInfo, registerInfo, h, ttrace] = ...
                                StackSet.LoadOneStimulus_ratio(ServerDir, Cam, p, ...
                                iStim, RepeatList, ResizeFactor, ratioInfo, ...
                                registerInfo, 'suffix', suffix);
                    else
                        S = StackSet.LoadOneStimulus(Cam, p, iStim, RepeatList, ResizeFactor, true);
                    end
                case 4 %21.01.14 AP
                    S = StackSet.LoadAllStimuliSelectedRepeat(Cam, p, RepeatList, ResizeFactor);
            end
            
            if isempty(S), return; end
            if newSave
                if ~strcmp(Cam.FileString , '_ratio')
                    S.SaveStacks(DIRS.Stacks, DIRS.Temp, Cam.FileString, suffix);
                
                elseif strcmp(Cam.FileString , '_ratio')%14/8/15
                    Sratio = S.GetOneCondition(1);
                    Ssum = S.GetOneCondition(2);
                
                    Sratio.SaveStacks(DIRS.Stacks, DIRS.Temp, '_ratio', suffix);
                    [StackDir, ThisFileName] = Ssum.SaveStacks(DIRS.Stacks, DIRS.Temp, '_sum', suffix);     
                    
                    %check if the sum stack is properly saved. added on
                    %2/9/2015
                    try
                        Ssum = tools.LoadMyStacks(StackDir, ThisFileName);
                    catch err
                        tools.SaveMyStacks(S, StackDir, [], ThisFileName);
                        disp(['Re-saving ' ThisFileName]);
                    end
                end
            end
            
            if exist('h','var') 
                StackDir = tools.getDirectory( ServerDir, p, ResizeFactor, iStim, RepeatList, Cam, suffix);
                tools.screen2pngDS(fullfile(StackDir, ['TtracePspec_' num2str(RepeatList)]), h);
                %save(fullfile(StackDir, ['TtracePspec_' num2str(RepeatList)]), 'ttrace');
                close(h);%15/5/14 DS
            end
            if exist('ttrace','var')
                StackDir = tools.getDirectory( ServerDir, p, ResizeFactor, iStim, RepeatList, Cam, suffix);
                save(fullfile(StackDir, ['TtracePspec_' num2str(RepeatList)]), 'ttrace');
            end
        otherwise
            S = StackSet;
            return
    end
end
