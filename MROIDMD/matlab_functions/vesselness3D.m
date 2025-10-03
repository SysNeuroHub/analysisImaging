function Vesselness = vesselness3D(V, sigmas, alpha, beta, c)

% V: 3D volume (double)
% sigmas: array of scales, e.g. [1 2 3 4]
% alpha, beta, c: Frangi parameters

V = double(V);
sz = size(V);

Vesselness = zeros(sz);

for s = sigmas
    %% 1. Create Gaussian derivative filters at scale s
    hsize = 2*ceil(3*s)+1;   % filter size
    [x,y,z] = ndgrid(-ceil(3*s):ceil(3*s));
    
    % Gaussian
    g = exp(-(x.^2 + y.^2 + z.^2)/(2*s^2));
    g = g / sum(g(:));
    
    % Derivatives of Gaussian (second order)
    Dxx = (x.^2/s^4 - 1/s^2) .* g;
    Dyy = (y.^2/s^4 - 1/s^2) .* g;
    Dzz = (z.^2/s^4 - 1/s^2) .* g;
    Dxy = (x.*y/s^4) .* g;
    Dxz = (x.*z/s^4) .* g;
    Dyz = (y.*z/s^4) .* g;
    
    %% 2. Convolve volume with filters to get Hessian components
    Hxx = convn(V, Dxx, 'same');
    Hyy = convn(V, Dyy, 'same');
    Hzz = convn(V, Dzz, 'same');
    Hxy = convn(V, Dxy, 'same');
    Hxz = convn(V, Dxz, 'same');
    Hyz = convn(V, Dyz, 'same');
    
    %% 3. Compute eigenvalues of Hessian
    Lambda1 = zeros(sz); Lambda2 = zeros(sz); Lambda3 = zeros(sz);
    
    % Reshape to vectorized form
    nvox = numel(V);
    H = [Hxx(:) Hxy(:) Hxz(:);
         Hxy(:) Hyy(:) Hyz(:);
         Hxz(:) Hyz(:) Hzz(:)];
    
    % Compute eigenvalues voxel by voxel (loop because eig doesn't vectorize easily)
    for i = 1:nvox
        Hi = [Hxx(i) Hxy(i) Hxz(i);
              Hxy(i) Hyy(i) Hyz(i);
              Hxz(i) Hyz(i) Hzz(i)];
        e = eig(Hi);
        % sort by absolute value
        [~,idx] = sort(abs(e));
        e = e(idx);
        Lambda1(i) = e(1);
        Lambda2(i) = e(2);
        Lambda3(i) = e(3);
    end
    
    Lambda1 = reshape(Lambda1, sz);
    Lambda2 = reshape(Lambda2, sz);
    Lambda3 = reshape(Lambda3, sz);
    
    %% 4. Vesselness measure
    Ra = abs(Lambda2) ./ abs(Lambda3);
    Rb = abs(Lambda1) ./ sqrt(abs(Lambda2 .* Lambda3));
    S  = sqrt(Lambda1.^2 + Lambda2.^2 + Lambda3.^2);
    
    Vn = exp(-(Ra.^2)/(2*alpha^2)) .* ...
         exp(-(Rb.^2)/(2*beta^2)) .* ...
         (1 - exp(-(S.^2)/(2*c^2)));
    
    % Keep only dark tubes
    Vn(Lambda2 > 0 | Lambda3 > 0) = 0;
    
    % Take maximum across scales
    Vesselness = max(Vesselness, Vn);
end
end
