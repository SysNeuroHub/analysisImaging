addpath('\\zserver\Code\stacks')

addpath '\\zserver.ioo.ucl.ac.uk\Code\Spikes'; % SetDefaultDirs, ProtocolLoad
SetDefaultDirs;
DIRS.Stacks = '\\zserver3.ioo.ucl.ac.uk\Data\Stacks';%where stackset data is saved

% input to Imager 
Magnification = 1.6;
HardwareBinning = 4;
FileString = '';

%% information on experimental condition
Exps.animal= 'M150501_SD';
Exps.iseries= 1;
Exps.iexp = 2;
Exps.ResizeFac= 0.5000;
Exps.Cam = tools.Imager('PCO', [], FileString, Magnification, HardwareBinning);
Exps.Cam.Info = 'ShutterMode:Rolling, ScanMode:slow, FrameRate:50';
Exps.Cam.DataDir = tools.getServerDirectory(Exps.animal);

iStim = 2;
iTr = 1;

%% information on stimulus protocol
p = ProtocolLoad(Exps);



%% one stimulus, single repeat
MyStack_oneRepeat = ...
    StackSet.LoadStacks( DIRS.Stacks, p, Exps.ResizeFac,iStim, ...
    iTr, Exps.Cam);

% avg repeats
saveDir = tools.getDirectory( DIRS.Stacks, p, Exps.ResizeFac, iStim, 1, Exps.Cam); 
MyStack_oneStim = tools.LoadMyStacks(saveDir,'nanMean');


% these are not integrated with timeline
% %% one stimulus, average repeats
% MyStack_oneStim2 = ...
%     StackSet.LoadStacks( DIRS.Stacks, p, Exps.ResizeFac,iStim, ...
%     [], Exps.Cam);
% 
% %% all stimuli, average repeats
% MyStack_allStim = ...
%     StackSet.LoadStacks( DIRS.Stacks, p, Exps.ResizeFac,[], ...
%     [], Exps.Cam);



