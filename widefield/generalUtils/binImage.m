function binnedImage = binImage(origImage, binFactor)
% binnedImage = binImage(origImage, binFactor)
% Bins image (e.g. 2x2 or 4x4) by taking the average of given pixels. 
% Inputs
%   origImage: image size Ypix x Xpix x nFrames. 
%   binFactor:B, for binning performed like BxB
% Outputs:
% binnedImage: image with nPix/B^2 pixels x nFrames. 
%
% Only works for integer B, probably.

% called in loadRawToDat.m

cFilt = ones(1, binFactor)/binFactor;
for im = 1:size(origImage,3)
    q = conv2(cFilt, cFilt, origImage(:,:,im), 'same');
    b = q(1:binFactor:end,1:binFactor:end);
    if im==1
        binnedImage = zeros(size(b,1), size(b,2), size(origImage,3), 'like', origImage);
    end
    binnedImage(:,:,im) = b;
end


% % Somehow this algorithm with convn is not faster, despite seeming more
% % natural to me. 
% cFilt = ones(binFactor, binFactor)/(binFactor^2);
% q = convn(origImage, cFilt, 'same');
% binnedImage = q(1:binFactor:end,1:binFactor:end,:);
