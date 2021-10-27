function regMovie = imageReg(movie, dx, dy, msize)
% imageReg(movie, dx, dy)
% returns a registered movie using the information of registration (dx, dy) (from RapidReg)
%
% imageReg(movie,dx,dy,msize)
% lets you specify margin size (0-1, DEFAULT:0)
%
% See also: tools.RapidReg, StackSet.Register

% 2014-1-29 DS made from rapidReg.m
% 2015-7-17 DS added window to image before applying fft

if nargin < 4
    msize = 0;
end
if isempty(msize)
    msize = 0;
end


disp('Making registered movie via imageReg...');

[h, w, nFrames] = size(movie);

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

%add mirror images to avoid fft artifact ripple - 2014.1.1 DS
orsize = size(movie);
xmsize = ceil(msize*orsize(2))-1; %margin size in x-axis
ymsize = ceil(msize*orsize(1))-1; %margin size in y-axis

lh = orsize(1) + 2*(ymsize+1);
lw = orsize(2) + 2*(xmsize+1);
lnr = ifftshift((-fix(lh/2):ceil(lh/2) - 1));
lnc = ifftshift((-fix(lw/2):ceil(lw/2) - 1));
[lnc, lnr] = meshgrid(lnc, lnr);


%%2D window to minimize boudary effect. 2015/7/17
w_func = @tukeywin;
wc=window(w_func,lh, 2*msize);
wr=window(w_func,lw, 2*msize);
[maskr,maskc]=meshgrid(wr,wc);

window2D = maskr.*maskc;


regMovie = zeros(h, w, nFrames, class(movie));
for t = 1:nFrames
    
    %     ftRegFrame = ftMovie(:,:,t).*exp(sqrt(-1)*2*pi*(-dy(t)*nr/h - dx(t)*nc/w));
    %     regMovie(:,:,t) = abs(ifft2(ftRegFrame));
    
    if msize>0
        lmovie = [flipdim(movie(:,1:xmsize+1,t),2) movie(:,:,t) flipdim(movie(:,orsize(2)-xmsize:end,t),2)];
        lmovie = [flipdim(lmovie(1:ymsize+1,:),1); lmovie; flipdim(lmovie(orsize(1)-ymsize:end,:),1)];
        lftMovie = fft2(lmovie.*window2D);
        
        ftRegFrame = lftMovie.*exp(sqrt(-1)*2*pi*(-dy(t)*lnr/lh - dx(t)*lnc/lw));
        lregMovie = abs(ifft2(ftRegFrame));
        regMovie(:,:,t) = lregMovie(ymsize+2:orsize(1)+ymsize+1,xmsize+2:orsize(2)+xmsize+1);
        
    elseif msize == 0
        ftMovie = fft2(movie(:,:,t));
        ftRegFrame = ftMovie.*exp(sqrt(-1)*2*pi*(-dy(t)*nr/h - dx(t)*nc/w));
        regMovie(:,:,t) = abs(ifft2(ftRegFrame));
    end
end


%% Convert the registered movie to its original type
regMovie = origTypeFun(regMovie);
