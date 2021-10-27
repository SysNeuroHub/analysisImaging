function S = LoadAllStimuliSelectedRepeat(Cam, p, iRepeat, ResizeFactor)
% Loads desired repeat for all stimuli from raw data
%
% MyStack = LoadAllStimuliSelectedRepeats(Cam.DataDir, p, RepeatList, ResizeFactor, Cam.Type, Cam.FileString)
% Input: - Cam, an Imager object
%        - p: p-file protocol;
%        - RepeatList: list of repeats you want to load and average
%        - ResizeFactor: resize the image stacks according to ResizeFactor
%
% Returns a StackSet with the values averaged across all repeats for all
% stimuli

if nargin<4 || isempty(ResizeFactor)   
    ResizeFactor = 1; 
end

fprintf('Loading averages across repeats for %s-%d-%d (%s, %d stimuli)\n',...
    p.animal, p.iseries, p.iexp, Cam.Type, p.nstim);
            
% load the first stimulus to initialize things
S = StackSet.LoadOneStimulus(Cam, p, 1, iRepeat, ResizeFactor );
if isempty(S), return; end

S.nConds = p.nstim;
S.Values = nan(S.nRows,S.nCols,S.nFrames,S.nConds);
            
% This loads one stimulus at a time, to save memory.
for iStim = 1:p.nstim
    StimAverage = zeros(S.nRows,S.nCols,S.nFrames); 
    nnnRepeats  = zeros(S.nRows,S.nCols,S.nFrames);
    ThisStack = StackSet.LoadOneStimulus(Cam, p,iStim, iRepeat, ResizeFactor );
    if ~isempty(ThisStack)
        StimAverage = nansum( cat(4,StimAverage, ThisStack.Values), 4 );
    end
    S.Values(:,:,:,iStim) = StimAverage;
end

S.Info.Stimulus = 'All';
S.Info.RepeatList = iRepeat;
S.Description   = sprintf('%s-%d-%d All Stimuli for selected repeat',p.animal,p.iseries,p.iexp);

end
