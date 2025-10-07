function niftiwrite_us(filename, scaleFactor)
%save a nii data after rescaling

if ~strcmp(filename(end-3:end),'.nii')
    error('niftiwrite_us only accepts nii file');
end

V = niftiread(filename);
info = niftiinfo(filename);

V = imresize3(V, scaleFactor, 'linear');
%info.Datatype = 'double';
info.PixelDimensions = info.PixelDimensions/scaleFactor;
info.Transform.T(1:3,1:3) = info.Transform.T(1:3,1:3)/scaleFactor;
info.ImageSize = size(V);

niftiwrite(V, [filename(1:end-4) '_us.nii'], info);
