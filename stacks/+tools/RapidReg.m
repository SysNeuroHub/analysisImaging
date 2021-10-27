function [dx, dy, target, regMovie, x, y] = RapidReg(movie, target, msize, precision, varargin)
% registers movie frames to a target frame using parallelisation
% and array vectorisation for efficiency
%
%   [DX, DY, TARGET, REG, X, Y] = RapidReg(MOVIE, TARGET) register
%   MOVIE (an (X,Y,T) array) to TARGET (either the target image frame (X,Y)
%   or 'auto', to find one automatically). Returns the registered movie, REG,
%   the x and y coordinates of the registered movies pixels relative to the
%   target (useful for when clipping is used), the translations required to
%   register each frame, DX and DY, and TAGRET, the target frame.
%   Optionally takes, 'noparallel', meaning use single-threaded codepath
%   instead of parallel, and 'clip', meaning clip the output image to keep
%   only pixels that have valid info at all times (i.e. haven't translated
%   outside target).
%
%   [...] = RapidReg(movie, target, msize)
%   allows to specify margin size to avoid edge artifact (0-1, DEFAULT:0)
%
%   [...] = RapidReg(movie, target, msize, precision)
%   allows to specify precision of the registration (DEFAULT:100). Unit of
%   the registration will be 1/precision [pixels].
%
%   [dx,dy,target] = RapidReg returns only the information needed to
%   register the images
%
%   [__] = RapidReg( __, Name, Value )
%   specifies additional information to register images using one or more
%   Name, Value pair arguments. Name, Value pair settings apply to all the
%   stacks loaded.
%   Name-Value pair arguments
%     'msize': size of the mirror image ([0 - 1])
%     'precision': 1/precision is minimum deviation estimated by dftregistration
%
% See also: tools.ImageReg, StackSet.Register

% 2013-07 CB created (heavily plagiarised from Mario Dipoppa's code)
% 2014-01 DS added varargin inputs

% TO DO: is msize necessary for RapidReg??


if nargout < 4
    GetInfoOnly = true;
else
    GetInfoOnly = false;
end
if nargin < 4
    precision = 100;
end
if isempty(precision)
    precision = 100;
end
if nargin < 3
    msize = 0;
end
if isempty(msize)
    msize = 0;
end

[h, w, nFrames] = size(movie);

fprintf('RapidReg..');

%% Setup
%convert movie data to an appropriate floating point type if necessary
dataType = class(movie);
switch dataType
    case {'int8' 'int16' 'uint8' 'uint16' 'int32' 'uint32'}
        %convert all integer types up to 32-bits to single
        movie = single(movie);
        origTypeFun = str2func(dataType);
    case {'int64'  'uint64'}
        %convert 64-bit integer types to double
        movie = double(movie);
        origTypeFun = str2func(dataType);
    case {'single' 'double'}
        %no conversion
        origTypeFun = str2func(dataType);%DS on 13.12.14
    otherwise
        error('''%s'' is not a recognised data type', dataType);
end

%create a Gaussian filter for filtering frames
hGauss = fspecial('gaussian', [5 5], 1);

%look for flag on whether to use parallel codepath
if any(cell2mat(strfind(varargin, 'nopar')) == 1)
    parallel = false;
else
    parallel = true;
end

if strcmpi(target, 'auto')
    %% Compute the best target frame
    fprintf('finding target..');
    %first compute a smoothed mean of each frame
    meanF = smooth(mean(reshape(movie, h*w, nFrames)));
    %now look in the middle third of the image frames for the minimum
    fromFrame = round(nFrames*1/3);
    toFrame = round(nFrames*2/3);
    [~, idx] = min(meanF(fromFrame:toFrame));
    minFrame = fromFrame + idx;
    %Gaussian filter the target image
    target = imfilter(movie(:,:,minFrame), hGauss, 'same', 'replicate');
    % AR added display output of which frame number is being used as target
    fprintf(['target frame is ',num2str(minFrame),'..']);
end

%% preprocessing target and movie
[ly, lx] = size(target);
maskSlope   = 2;%0.5;%1.2; % 9/8/16
smoothSigma = 2; %1.15 %9/8/16

% smoothing filter in frequency domain ... only for target image
hgx = exp(-(((0:lx-1) - fix(lx/2))/smoothSigma).^2);
hgy = exp(-(((0:ly-1) - fix(ly/2))/smoothSigma).^2);
hg = hgy'*hgx;
fhg = real(fftn(ifftshift(single(hg/sum(hg(:))))));


% Taper mask
[ys, xs] = ndgrid(1:ly, 1:lx);
ys = abs(ys - mean(ys(:)));
xs = abs(xs - mean(xs(:)));
mY      = max(ys(:)) - 4;%round(0.1*ly);
mX      = max(xs(:)) - 4;%round(0.1*lx);
maskMul = double(1./(1 + exp((ys - mY)/maskSlope)) ./(1 + exp((xs - mX)/maskSlope)));
maskOffset = mean(target(:))*(1 - maskMul);

movie = bsxfun(@plus, maskOffset, bsxfun(@times, maskMul, movie));%heavy??


%% Fourier transform the movie frames, unfiltered and filtered
fprintf('filtering..');
ftMovie = fft2(movie);
fttarget = fft2(target);
%fttarget = conj(fftn(target)); %NG

%phase correlation
eps0 = double(1e-20);
fttarget = fttarget ./ (abs(fttarget) + eps0).*fhg;
ftMovie = ftMovie ./ (abs(ftMovie) + eps0);


%% Compute required displacement and register each frame
dx = zeros(1, nFrames);
dy = zeros(1, nFrames);
nr = ifftshift((-fix(h/2):ceil(h/2) - 1));
nc = ifftshift((-fix(w/2):ceil(w/2) - 1));
[nc, nr] = meshgrid(nc, nr);
regMovie = zeros(h, w, nFrames, class(movie));
fprintf('registering..');

if parallel
    %% Register in parallel
    temporaryPool = isempty(gcp('nocreate'));
    if temporaryPool
        parpool;%create default worker pool
    end
    try
        %do parallel loops in chunks of data to prevent matlab choking
        chunkSize = min(14000, nFrames); %frames
        nChunks = ceil(nFrames/chunkSize);
        for i = 0:(nChunks - 1)
            sidx = i*chunkSize + 1;
            eidx = min((i + 1)*chunkSize, nFrames);
            parfor t = sidx:eidx
                %find the best registration translation
                output = dftregistration(fttarget, ftMovie(:,:,t), precision);
                dx(t) = output(4);
                dy(t) = output(3);
                
                if ~GetInfoOnly
                    %translate the original (i.e. unfiltered) frame
                    ftRegFrame = ftMovie(:,:,t).*exp(sqrt(-1)*2*pi*(-dy(t)*nr/h - dx(t)*nc/w));
                    regMovie(:,:,t) = abs(ifft2(ftRegFrame));
                end
            end
        end
        if temporaryPool
            delete(gcp('nocreate')); %close worker pool
        end
    catch ex
        if temporaryPool
            %in case of error, ensure temporary worker pool is closed
            delete(gcp('nocreate'));
        end
        rethrow(ex)
    end
else
    %% Register sequentially
    
    %add mirror images to avoid fft artifact ripple - 2014.1.1 DS
    orsize = size(target);
    xmsize = ceil(msize*orsize(2))-1; %margin size in x-axis
    ymsize = ceil(msize*orsize(1))-1; %margin size in y-axis
    
    
    
    %[lh, lw, ~] = size(lmovie);
    lh = orsize(1) + 2*(ymsize+1);
    lw = orsize(2) + 2*(xmsize+1);
    lnr = ifftshift((-fix(lh/2):ceil(lh/2) - 1));
    lnc = ifftshift((-fix(lw/2):ceil(lw/2) - 1));
    [lnc, lnr] = meshgrid(lnc, lnr);
    
    
    for t = 1:nFrames
        %find the best registration translation
        output = dftregistration(fttarget, ftMovie(:,:,t), precision);
        dx(t) = output(4);
        dy(t) = output(3);
        %translate the original (i.e. unfiltered) frame
        
        if ~GetInfoOnly
            if msize>0
                lmovie = [flipdim(movie(:,1:xmsize+1,t),2) movie(:,:,t) flipdim(movie(:,orsize(2)-xmsize:end,t),2)];
                lmovie = [flipdim(lmovie(1:ymsize+1,:),1); lmovie; flipdim(lmovie(orsize(1)-ymsize:end,:),1)];
                lftMovie = fft2(lmovie);
                
                ftRegFrame = lftMovie.*exp(sqrt(-1)*2*pi*(-dy(t)*lnr/lh - dx(t)*lnc/lw));
                lregMovie = abs(ifft2(ftRegFrame));
                regMovie(:,:,t) = lregMovie(ymsize+2:orsize(1)+ymsize+1,xmsize+2:orsize(2)+xmsize+1);
                
            elseif msize == 0
                ftMovie = fft2(movie(:,:,t));
                ftRegFrame = ftMovie.*exp(sqrt(-1)*2*pi*(-dy(t)*nr/h - dx(t)*nc/w));
                regMovie(:,:,t) = abs(ifft2(ftRegFrame));
            end
        end
        
    end
end

%% If requested, clip the frames to the maximum fully valid region
if any(cell2mat(strfind(varargin, 'clip')) == 1)
    fprintf('clipping..');
    dxMax = max(0, ceil(max(dx)));
    dxMin = min(0, floor(min(dx)));
    dyMax = max(0, ceil(max(dy)));
    dyMin = min(0, floor(min(dy)));
    x = (1 + dxMax):(h + dxMin);
    y = (1 + dyMax):(h + dyMin);
    regMovie = regMovie(y,x,:);
else
    x = 1:w;
    y = 1:h;
end
%% Convert the registered movie to its original type
regMovie = origTypeFun(regMovie);

fprintf('.done\n');

end

