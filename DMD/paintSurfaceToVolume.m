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

if (volSize(1) ~= size(surfData,1)) || (volSize(2) ~= size(surfData,2))
    error('volume size does not match image size');
end

TexVol = zeros(volSize, 'like', surfData);


%surfData = flipud(surfData);
[rows, cols] = size(surfDepth);

% convert subscripts to linear indices
iy = (1:rows)' * ones(1,cols);
ix = ones(rows,1) * (1:cols);
%linIdx = sub2ind(volSize, ix(:), round(surfDepth(:)), iy(:));
linIdx = sub2ind(volSize, iy(:), ix(:), round(surfDepth(:)));

% flatten everything
linIdx = linIdx(:);
vals   = surfData(:);

% remove invalid (zero or out-of-bounds) indices
valid = linIdx > 0 & linIdx <= numel(TexVol) & surfDepth(:) > 0;
TexVol(linIdx(valid)) = vals(valid);
%TexVol = flip(TexVol,3); %yes i need this
end


%% another implementation 
% 
% 
% [X,Y] = meshgrid(1:size(surfDepth,2),1:size(surfDepth,1));
% Z = surfDepth;
% C = surfData;
% 
% h = surf(X, Y, Z, C);
% 
% X = h.XData;
% Y = h.YData;
% Z = h.ZData;
% C = h.CData;
% 
% % xv = linspace(min(X(:)), max(X(:)), Nx);   % Nx = desired number of voxels in x
% % yv = linspace(min(Y(:)), max(Y(:)), Ny);
% % zv = linspace(min(Z(:)), max(Z(:)), Nz);
% 
% xv = min(X(:)):max(X(:));
% yv = min(Y(:)):max(Y(:));
% zv = min(Z(:)):max(Z(:));
% 
% [Xv, Yv, Zv] = ndgrid(xv, yv, zv);
% 
% % option0: too slow
% % F = scatteredInterpolant(X(:), Y(:), Z(:), C(:), 'nearest', 'none');
% % V = F(Xv, Yv, Zv);
% 
% % option1: too slow
% % V1 = griddata(X, Y, Z, C, Xv, Yv, Zv, 'nearest');
% 
% 
% % option2
% Nx = volSize(1);
% Ny = volSize(3);
% Nz = volSize(2);
% TexVol = zeros(Nx, Ny, Nz);
% ix = round( rescale(X, 1, Nx) );
% iy = round( rescale(Y, 1, Ny) );
% iz = round( rescale(Z, 1, Nz) );
% 
% linIdx = sub2ind(size(TexVol), ix, iy, iz);
% TexVol(linIdx) = C;  % assign color at those voxels
% 
% TexVol = permute(TexVol, [1 3 2]);