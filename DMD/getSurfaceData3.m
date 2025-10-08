function [surfDepth, surfData] = getSurfaceData3(mr_brain, direction, threshold, smoothSigma)
% getSurfaceData  Extracts a smooth surface map and intensity values from a 3D volume,
%                 using only voxel-based operations (no isosurface).
%
%   [surfDepth, surfData] = getSurfaceData(mr_brain)
%   [surfDepth, surfData] = getSurfaceData(mr_brain, direction, threshold, smoothSigma)
%
%   Inputs:
%       mr_brain   - 3D matrix (e.g. MRI volume)
%       direction  - 'first' (top surface) or 'last' (bottom surface) [default: 'first']
%       threshold  - intensity threshold to define object [default: 0.2 * max(mr_brain(:))]
%       smoothSigma - Gaussian smoothing of surface depth [default: 2]
%
%   Outputs:
%       surfDepth  - depth index (z) of detected surface per (x,y)
%       surfData   - voxel intensity at that surface

% --- Defaults ---
if nargin < 2 || isempty(direction), direction = 'first'; end
if nargin < 3 || isempty(threshold), threshold = 0.2 * max(mr_brain(:)); end
if nargin < 4 || isempty(smoothSigma), smoothSigma = []; end

% --- Consistent orientation (match your original) ---
aa_brain = permute(mr_brain,[1 3 2]);  
[rows, cols, slices] = size(aa_brain);

% --- Binary mask of object ---
mask = aa_brain > threshold;

% --- Find surface depth voxelwise ---
switch direction
    case 'first'  % topmost voxel (smallest z)
        mask_cumsum = cumsum(mask,3);
        surfaceMask = (mask_cumsum == 1) & mask;
        [~, surfZ] = max(surfaceMask, [], 3);

    case 'last'   % deepest voxel (largest z)
        mask_flip = flip(mask,3);
        mask_cumsum = cumsum(mask_flip,3);
        surfaceMask = (mask_cumsum == 1) & mask_flip;
        [~, surfZ] = max(surfaceMask, [], 3);
        surfZ = slices - surfZ + 1; % invert z index

    otherwise
        error('direction must be ''first'' or ''last''.');
end

if ~isempty(smoothSigma)
    % --- Smooth the surface depth map ---
    surfZ = imgaussfilt(surfZ, smoothSigma);
end

% --- Sample data on the detected surface ---
[xGrid, yGrid] = ndgrid(1:rows, 1:cols);

% --- Rotate back for visualization consistency ---
surfDepth = double(fliplr(rot90(surfZ)));


if nargout>1
    surfD = interp3(aa_brain, yGrid, xGrid, surfZ, 'linear', 0);
    surfData  = double(fliplr(rot90(surfD)));
end
end
