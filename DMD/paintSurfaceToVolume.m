function TexVol = paintSurfaceToVolume(surfDepth, surfData, volSize)
% paintSurfaceToVolume  Create 3D volume assigning surface data to detected surface voxels.
%
%   TexVol = paintSurfaceToVolume(surfDepth, surfData, volSize)
%
%   Inputs:
%       surfDepth - matrix of z indices (e.g., from getSurfaceData)
%       surfData  - matrix of voxel intensities at surface
%       volSize   - 3-element vector, e.g. size(mr_brain)
%
%   Output:
%       TexVol    - 3D volume with surfData placed at [x, y, z = surfDepth(x,y)]

if (volSize(1) ~= size(surfData,2)) || (volSize(3) ~= size(surfData,1))
    error('volume size does not match image size');
end

TexVol = zeros(volSize, 'like', surfData);


surfData = flipud(surfData);
[rows, cols] = size(surfDepth);

% convert subscripts to linear indices
iy = (1:rows)' * ones(1,cols);
ix = ones(rows,1) * (1:cols);
linIdx = sub2ind(volSize, ix(:), round(surfDepth(:)), iy(:));

% flatten everything
linIdx = linIdx(:);
vals   = surfData(:);

% remove invalid (zero or out-of-bounds) indices
valid = linIdx > 0 & linIdx <= numel(TexVol) & surfDepth(:) > 0;
TexVol(linIdx(valid)) = vals(valid);
%TexVol = flip(TexVol,3);
end
