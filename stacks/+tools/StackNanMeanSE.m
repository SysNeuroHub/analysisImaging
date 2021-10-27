function [S_nanMean, S_nanSE, validFrames] = StackNanMeanSE(ServerDir, Cam, p, iStim, RepeatList, ResizeFactor, suffix)
% [S_nanMean, S_nanSE] = StackNanMeanSE(...) returns StackSet mean and SE across
% repeats, skipping frames that were interpolated with fixMissingFramesFunc.
%
% [S_nanMean, S_nanSE, validFrames] = StackNanMeanSE(...) also returns
% frame indexes that were used for each repeat
%
% Currently only applicable to stacks of LoadOption 3 (single stimulus, specified repeats)
%
% iStim: scalar value
% RepeatList: vector or empty to specify repeat index to process
% nanMean/nanSE

% TODO: make a summary figure (timecourse, powerspectrum, amp/phase map) 
% change name of this function? since it is not really nanmean, or add
% option to calculate nanmean?

% 2014-12-16 DS added exception when interp information is not available
% 2015-5-1 DS added 7th input (suffix)

% this function will be replaced by StackEvRepeats in future

if nargin<7
    suffix = '';
end

if isempty(RepeatList)
    RepeatList = 1:p.nrepeats;
end

SumStack = [];
SquareStack = [];
validFrames = [];

for iTr = RepeatList
    MyStack = StackSet.LoadStacks( ServerDir, p, ResizeFactor,iStim, iTr, Cam, 'suffix',suffix);%10/7/15   
    %load info
    [StackDir, ThisFileName] = tools.getDirectory( ServerDir, p, ResizeFactor, iStim, iTr, Cam, suffix);
    
    if exist([StackDir '/info' ThisFileName], 'file')
        
        information = load([StackDir '/info' ThisFileName]);
        information = information.info;
        
        if strcmp(Cam.FileString, '_ratio')
            
            if isfield(information, 'ratioInfo') && isfield(information.ratioInfo, 'interp')
                missedFrames = information.ratioInfo.interp;
            else
                missedFrames = [];
            end
        else %cam1 or cam2
            if isfield(information, 'timeInfo') && isfield(information.timeInfo, 'interp')
                missedFrames = information.timeInfo.interp;
            else
                missedFrames = [];
            end
        end
    else
        display('Information for missed frames was not found.');
        missedFrames = [];
    end
    
    
    if ~exist('TimeStamp_align', 'var')
        TimeStamp_align = MyStack.TimeVec;
    end
        
    [MyStack.Values, idx_align] = tools.alignTimeStamp(MyStack.Values, MyStack.TimeVec, TimeStamp_align);
    [~, missedFrames] = intersect(idx_align, missedFrames);
    MyStack.nFrames = size(MyStack.Values,3);
    
    nmissedFrames = setxor(1:size(MyStack.Values,3), missedFrames);
    if ~isempty(SumStack)
        nmissedFrames = intersect(1:size(SumStack, 3), nmissedFrames);
    end
%     nmissedTime = ones(size(MyStack.Values,3), 1);
    nmissedTime = ones(1, size(MyStack.Values,3)); %26/4/14
    nmissedTime(missedFrames) = 0;
    
    MyStack.Values(:,:,missedFrames) = 0;
    
    if isempty(SumStack)
        SumStack = zeros(size(MyStack.Values), 'single'); %16/5/14 double > single
    end
    if isempty(SquareStack)
        SquareStack = zeros(size(MyStack.Values), 'single');  %16/5/14 double > single
    end
    if isempty(validFrames)
        validFrames = zeros(p.nrepeats, size(MyStack.Values, 3));
    end 
    tsize = min(size(SumStack,3),size(MyStack.Values,3));
    SumStack(:,:,nmissedFrames) = SumStack(:,:,nmissedFrames) + MyStack.Values(:,:,nmissedFrames);
    SquareStack(:,:,nmissedFrames) = SquareStack(:,:,nmissedFrames) + MyStack.Values(:,:,nmissedFrames).^2;
    validFrames(iTr, 1:min(size(MyStack.Values,3), length(nmissedTime))) ...
        = nmissedTime(1:min(size(MyStack.Values,3), length(nmissedTime)));
    
end
FrameRate = MyStack.FrameRate;

MeanStack = zeros(size(SumStack,1),size(SumStack,2),size(SumStack,3),'single'); %16/5/14
SEStack = zeros(size(SumStack,1),size(SumStack,2),size(SumStack,3),'single'); %16/5/14
for tt = 1:size(SumStack,3)
    numTr = sum(validFrames(:,tt));
    MeanStack(:,:,tt) = SumStack(:,:,tt) / numTr;
    SEStack(:,:,tt) = sqrt(SquareStack(:,:,tt)/numTr - MeanStack(:,:,tt).^2) / sqrt(numTr);
end
clear SumStack;

MyStack.Values = [];
S = MyStack;
S.nConds = 1;
S.TimeVec = TimeStamp_align;
S.nFrames = length(TimeStamp_align);

S.Values = MeanStack;
S.Info.RepeatList = 'nanMean';
S.Description = sprintf('nan mean across %d repeats [%s] ',  length(RepeatList),num2str(RepeatList));
S_nanMean = S;

S.Values = SEStack;
S.Info.RepeatList = 'nanSE';
S.Description = sprintf('nan SE across %d repeats [%s] ',  length(RepeatList),num2str(RepeatList));
S_nanSE = S;



