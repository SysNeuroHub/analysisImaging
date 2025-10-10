function [surfX, surfY, surfZ,  S] = vol2Surf(V, th_depth)

% Extract the surface voxels
% Get surface voxels using isosurface
fv = isosurface(V, 0.5);

% Compute surface normals
normals = isonormals(V, fv.vertices);
normals = normals ./ vecnorm(normals, 2, 2);

% Keep surfaces facing +z direction and depth>50
%keepIdx = find((normals(:,1) > 0).*(fv.vertices(:,1) > th_depth));

% % depths = unique(fv.vertices(:,1));
% % for dd = 1:numel(depths)
% % fv.vertices(:,1) == depths()

fv.vertices = fv.vertices(keepIdx,:);

% Round coordinates to voxel indices
surfX = round(fv.vertices(:,1)); %round(fv.vertices(:,2));% dorso-ventral
surfY = round(fv.vertices(:,2)); %round(fv.vertices(:,3));% left-right
surfZ = round(fv.vertices(:,3)); %round(fv.vertices(:,1));% antero-posterior

if nargout > 3
    S = zeros(size(V));
    idx = sub2ind(size(V), surfY, surfX, surfZ);
    S(idx) = 1;
end