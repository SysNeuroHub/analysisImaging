function S = LoadAllStimuliAverageRepeats(Cam, p, ResizeFactor)
% Average all repeats to all stimuli from raw data
%
% MyStack = LoadAllStimuliAverageRepeats(Cam.DataDir, p, ResizeFactor, Cam.Type, Cam.FileString)
% Input: - Cam, an Imager object
%        - p: p-file protocol;
%        - ResizeFactor: resize the image stacks according to ResizeFactor
%
% Returns a StackSet with the values averaged across all repeats for all
% stimuli
%
% See also: tools.Imager, StackSet.LoadOneStimulus

if nargin<3 || isempty(ResizeFactor)   
    ResizeFactor = 1; 
end

fprintf('Loading averages across repeats for %s-%d-%d (%s, %d stimuli)\n',...
    p.animal, p.iseries, p.iexp, Cam.Type, p.nstim);
            
% load the first stimulus, first repeat, to initialize things
S = StackSet.LoadOneStimulus(Cam, p, 1, 1, ResizeFactor );
if isempty(S), return; end

S.nConds = p.nstim;
S.Values = nan(S.nRows,S.nCols,S.nFrames,S.nConds);
            
% This loads one stimulus, one repeat at a time, to save memory.
for iStim = 1:p.nstim
    StimAverage = zeros(S.nRows,S.nCols,S.nFrames); 
    nnnRepeats  = zeros(S.nRows,S.nCols,S.nFrames);
    for iRepeat = 1:p.nrepeats
        ThisStack = StackSet.LoadOneStimulus(Cam, p,iStim, iRepeat, ResizeFactor );
        if ~isempty(ThisStack)
            StimAverage = nansum( cat(4,StimAverage, ThisStack.Values), 4 );
            nnnRepeats  = nnnRepeats + isfinite( ThisStack.Values );
        end
    end
    S.Values(:,:,:,iStim) = StimAverage ./ nnnRepeats;
end

S.Info.Stimulus = 'All';
S.Info.RepeatList = 'Average';
S.Description   = sprintf('%s-%d-%d All Stimuli Average across repeats',p.animal,p.iseries,p.iexp);

end
