% ************************************************************************
% average2Filter
%   2D average filter filter
%   filter uses zero padding and convolution
%   INPUTs:
%   inputSequence: sequence of data
%   window: size of averaging kernel
%   passes
%   firstImage: first image to filter
%   lastImage: last image to filter
%
% WALTHER, March 2010
% 2014-02-13 DS changed nargin order. sigma 7th > 4th, type 6th > 3rd
% remove firstFrame, lastFrame
% 2017-11-05 added type = 'median'
% WILL BE REMOVED IN FUTURE
% ************************************************************************

function sequence = averageFilter2(inputSequence, window, type, sigma, passes, imgtlbx)

switch nargin
    case 0
        error('input missing');
    case 1
        window = 5;
        passes = 1;
        type = 'average';
        imgtlbx = true;
        sigma = 1.5;
    case 2
        passes = 1;
        type = 'average';
        imgtlbx = true;
        sigma = 1.5;
    case 3
        passes = 1;
        imgtlbx = true;
        sigma = 1.5;
    case 4
        passes = 1;
        imgtlbx = true;
    case 5
        imgtlbx = true;
end

frames = size(inputSequence,3);
firstFrame = 1;
lastFrame = size(inputSequence, 3);

window = max(1, window);
window = min(window, 50);

passes = max(1, passes);
passes = min(passes, 20);

limit = [firstFrame; lastFrame];
limit( limit > frames ) = frames;
limit( limit < 1 ) = 1;
limit = sort( limit );
firstFrame = limit(1);
lastFrame = limit(2);
clear limit;

if isinteger(inputSequence)
    inputSequence = single(inputSequence);
end
sequence = single(inputSequence);

if strcmp(type, 'average')
    kernel = fspecial(type, window);
elseif strcmp(type, 'gaussian')
    kernel = fspecial(type, window, sigma);
end

if strcmp(type, 'average') || strcmp(type, 'gaussian')
    
    if imgtlbx
        for k = firstFrame:lastFrame
            for repeat = 1:passes
                sequence(:,:,k) = single(imfilter(inputSequence(:,:,k), kernel, 'symmetric'));
            end
        end
    else
        for k = firstFrame:lastFrame
            for repeat = 1:passes
                sequence(:,:,k) = single(filter2(kernel, inputSequence(:,:,k), 'same'));
            end
        end
    end
    
elseif strcmp(type, 'median')
    for k = firstFrame:lastFrame
        sequence(:,:,k) = single(medfilt2(inputSequence(:,:,k), [window window], 'symmetric'));
    end
end