
function sequence = differentialSequence(inputSequence, firstBaseFrame, lastBaseFrame, firstFrame, lastFrame, meanImage)
%   differentialSequence(inputSequence, firstBaseFrame, lastBaseFrame)
%   produces a differential image sequence
% 
%   INPUTs:
%   inputSequence: matlab sequence of images
%   firstBaseFrame: first frame of baseline
%   lastBaseLine: last frame of baseline
% 
%   See also: StackSet.differential, StackSet.LoadOneStimulus_ratio

%   Walther, March 2010


switch nargin
    case {0, 1, 2}
        error('input error!')
    case 3
        firstFrame = 1;
        lastFrame = size(inputSequence, 3);
    case 4
        lastFrame = size(inputSequence, 3);
end

[height, width, frames] = size(inputSequence);
sequence = ones(height, width, frames, 'single');%ok??

base = [ firstBaseFrame; lastBaseFrame ];
base( base > frames ) = frames;
base( base < 1 ) = 1;
base = sort( base);

limit = [firstFrame; lastFrame];
limit( limit > frames ) = frames;
limit( limit < 1 ) = 1;
limit = sort( limit );

% baseline image
if nargin < 6 || isempty(meanImage)
    meanImage = single(mean(inputSequence(:,:,base(1):base(2)), 3));
end

% differential sequence

for k=limit(1):limit(2)
    sequence(:,:,k) = single(inputSequence(:,:,k))./ meanImage;
end

end


