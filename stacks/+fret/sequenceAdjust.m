function sequenceCorrected = sequenceAdjust( sequence, correction, firstFrameBase, lastFrameBase, firstFrame, lastFrame, sequenceBaseImage)
% returns sequence adjusted by noise-gain equalization factor, for each pixel in space
%
% sequenceCorrected = sequenceAdjust( sequence, correction, firstFrameBase, lastFrameBase, firstFrame, lastFrame)
% INPUTS:
%   sequence: x-y-t tensor
%   correction: noise-gain equalization factor, either scalar or 2D matrix
%   first(last)FrameBase: period to calculate the base image, to which
%   correction factor is applied
%   first(last)Frame: period to calculate the noise-gain equalization
%
% See also: fret.makeratio, fret.slavemasterAdjustFactor

% 2014-9-4 DS renamed sequenceAdjust2D to sequenceAdjust, 
% 2015-7-6 DS added 7th input sequenceBaseImage (3rd & 4th inputs are ignored)

switch nargin
    case {0, 1}
        error('SEQUENCE ADJUST: input missing');
    case 2
        firstFrameBase = 1;
        lastFrameBase = inf;
        firstFrame = 1;
        lastFrame = inf;
    case 3
        lastFrameBase = inf;
        firstFrame = 1;
        lastFrame = inf;
    case 4
        firstFrame = 1;
        lastFrame = inf;
    case 5
        lastFrame = inf;
end

frames = size(sequence,3);

baseFrame = [firstFrameBase; lastFrameBase];
baseFrame( baseFrame > frames ) = frames;
baseFrame( baseFrame < 1 ) = 1;
baseFrame = sort( baseFrame );

limitFrame = [firstFrame; lastFrame];
limitFrame( limitFrame > frames ) = frames;
limitFrame( limitFrame < 1 ) = 1;
limitFrame = sort( limitFrame );

% adjust sequence

sequenceCorrected = sequence;
if nargin < 7
    sequenceBaseImage = single(mean(sequence(:,:,baseFrame(1):baseFrame(2)), 3));
end

for k = limitFrame(1):limitFrame(2)
    %just modified to ".*" on 6.22.2012
      sequenceCorrected(:,:,k) = uint16(correction .* (single(sequence(:,:,k)) - sequenceBaseImage) + sequenceBaseImage);
end

end

