function ratio = ratioOfSequences(inputSequence, referenceSequence , firstFrame, lastFrame)
% ratio = ratioOfSequences(inputSequence, referenceSequence , firstFrame, lastFrame)
% returns the ratio of referenceSequence divided by inputSequence
%
% Inputs:
%   firstFrame as Integer: first sequence to be averaged
%   lastFrame as Integer: last sequence to be averaged
%
% Outputs:
%   average as Single: averaged image
%
% See also: fret.makeratio
%
% WALTHER March 2010
% 2014 DS converted output to Single precision


[height, width, frames] = size(referenceSequence);

switch nargin
    case {0, 1}
        error('input missing');
    case 2
        firstFrame = 1;
        lastFrame = inf;
    case 3
        lastFrame = inf;
end

if (size(inputSequence,1) ~= height) || (size(inputSequence, 2) ~= width) || (size(inputSequence, 3) ~= frames)
    error('input/reference mismatch');
end

limit = [firstFrame; lastFrame];
limit( limit > frames ) = frames;
limit( limit < 1 ) = 1;
limit = sort( limit );

ratio = ones(height, width, frames, 'single');

for k = limit(1):limit(2)
%     if registration
%             inputSequence(:,:,k) = imtransform(inputSequence(:,:,k), imageTransform,...
%                 'FillValues', inf, 'XData', [1, width], 'YData', [1, height]); 
%     end
    ratio(:,:,k) = single(referenceSequence(:,:,k)) ./ single(inputSequence(:,:,k));
end

end

