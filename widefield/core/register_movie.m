function dreg = register_movie(data, ops, ds)
% dreg = register_movie(data, ops, ds)
% returns a movie after image registration according to ds
% in this function, data size is not changed 

% called in :
% align_iterative for determineTargetFrame
% registerDatFile to save registered movie to dat file

if ops.useGPU
    data = gpuArray(single(data));
end
[Ly, Lx, NT] = size(data);

Ny = ifftshift([-fix(Ly/2):ceil(Ly/2)-1]);
Nx = ifftshift([-fix(Lx/2):ceil(Lx/2)-1]);
[Nx,Ny] = meshgrid(Nx,Ny);
Nx = Nx / Lx;
Ny = Ny / Ly;


if strcmp(ops.RegPrecision, 'same')
   dreg  = zeros(size(data), 'like', data);
else
    if ops.useGPU
        dreg = gpuArray.zeros(size(data), ops.RegPrecision);
    else
        dreg = zeros(size(data), ops.RegPrecision);
    end
end

if ops.useGPU
    ds = gpuArray(ds);
    Nx = gpuArray(single(Nx));
    Ny = gpuArray(single(Ny));
end

for i = 1:NT
    dph         = 2*pi*(ds(i,1)*Ny + ds(i,2)*Nx);
    fdata       = fft2(single(data(:,:,i)));
    dreg(:,:,i) = real(ifft2(fdata .* exp(1i * dph)));
end

if ops.useGPU
    dreg = gather(dreg);
end
