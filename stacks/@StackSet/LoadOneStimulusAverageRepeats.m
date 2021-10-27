function S = LoadOneStimulusAverageRepeats(Cam, p, ResizeFactor, iStim)
% Average all repeats to one stimulus from raw data
%
% MyStack = LoadOneStimulusAverageRepeats(Cam.DataDir, p, ResizeFactor, Cam.Type, Cam.FileString)
% Input: - Cam, an Imager object
%        - p: p-file protocol;
%        - ResizeFactor: resize the image stacks according to ResizeFactor
%
% Returns a StackSet with the values averaged across all repeats for single
% stimulus
% 
% See also: tools.Imager, StackSet.LoadOneStimulus


if nargin<3 || isempty(ResizeFactor)
    ResizeFactor = 1;
end

if length(iStim) > 1
    error('4th input(iStim) needs to be a scalar');
end

fprintf('Loading the first stimulus, first repeat, to initialize things\n');%need to be changed
S = StackSet.LoadOneStimulus(Cam, p, iStim, 1, ResizeFactor );
if isempty(S), return; end

S.nConds = 1;
S.Values = [];
% S.Values = nan(S.nRows,S.nCols,S.nFrames);

fprintf('Loading averages across repeats for %s-%d-%d (%s, stimuli)\n',...
    p.animal, p.iseries, p.iexp, Cam.Type);

% This loads one stimulus, one repeat at a time, to save memory.
StimAverage = zeros(S.nRows,S.nCols,S.nFrames);
if ResizeFactor == 1
    StimAverage = uint32(StimAverage);
end
nnnRepeats  = uint8(zeros(S.nRows,S.nCols,S.nFrames));
for iRepeat = 1:p.nrepeats
    ThisStack = StackSet.LoadOneStimulus(Cam, p,iStim, iRepeat, ResizeFactor );%return double due to imresizing
    
    
    if ResizeFactor == 1;
        ThisStack.Values = int32(ThisStack.Values);
    end
    if ~isempty(ThisStack)
        if ResizeFactor == 1;
            StimAverage = int32(nansum( cat(4,StimAverage, ThisStack.Values), 4 ));
        else
            StimAverage = nansum( cat(4,StimAverage, ThisStack.Values), 4 );
        end
        nnnRepeats  = nnnRepeats + uint8(isfinite( ThisStack.Values ));
    end
end
S.Values = double(StimAverage) ./ double(nnnRepeats);


S.Info.Stimulus = iStim;
S.Info.RepeatList = 'Average';
S.Description   = sprintf('%s-%d-%d Stim %d Average across repeats',...
    p.animal,p.iseries,p.iexp, iStim);

end
