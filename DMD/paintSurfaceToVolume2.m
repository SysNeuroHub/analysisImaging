function TexVol = paintSurfaceToVolume2(surfDepth, surfData, volSize)

[X,Y] = meshgrid(1:size(surfDepth,2),1:size(surfDepth,1));
Z = surfDepth;
C = surfData;

h = surf(X, Y, Z, C);

X = h.XData;
Y = h.YData;
Z = h.ZData;
C = h.CData;

% xv = linspace(min(X(:)), max(X(:)), Nx);   % Nx = desired number of voxels in x
% yv = linspace(min(Y(:)), max(Y(:)), Ny);
% zv = linspace(min(Z(:)), max(Z(:)), Nz);

xv = min(X(:)):max(X(:));
yv = min(Y(:)):max(Y(:));
zv = min(Z(:)):max(Z(:));

[Xv, Yv, Zv] = ndgrid(xv, yv, zv);

% option0: too slow
% F = scatteredInterpolant(X(:), Y(:), Z(:), C(:), 'nearest', 'none');
% V = F(Xv, Yv, Zv);

% option1: too slow
% V1 = griddata(X, Y, Z, C, Xv, Yv, Zv, 'nearest');


% option2
Nx = volSize(1);
Ny = volSize(3);
Nz = volSize(2);
TexVol = zeros(Nx, Ny, Nz);
ix = round( rescale(X, 1, Nx) );
iy = round( rescale(Y, 1, Ny) );
iz = round( rescale(Z, 1, Nz) );

linIdx = sub2ind(size(TexVol), ix, iy, iz);
TexVol(linIdx) = C;  % assign color at those voxels

TexVol = permute(TexVol, [1 3 2]);