function frameRecon = svdFrameReconstruct(U, V)
% frameRecon = svdFrameReconstruct(U, V)
% returns 3D tensor (x,y,frames) reconstructed from U and V
% 
% Inputs:
% U: Y x X x nSVD
% V: nSVD x nFrames
% 
% Outputs:
% frameRecon: Y x X x nframes

% reshape U to be nPix x nSVD
Ur = reshape(U, size(U,1)*size(U,2), size(U,3));

% multiply and reshape back into Y x X
frameRecon = reshape(Ur*V, size(U,1), size(U,2), size(V,2));