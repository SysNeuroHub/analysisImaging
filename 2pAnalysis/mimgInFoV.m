function mimgOriginalSize = mimgInFoV(dat)
%mimgOriginalSize = mimgInFoV(dat)
% returns mean image (mimg) in FoV during recording
%16/7/20 created from Maruis2Craig

%Ly_temp = size(dat.mimg,1); %NG ops.mimg
%Lx_temp = size(dat.mimg,2); %NG ops.mimg
%mimg_c = zeros(Ly_temp,Lx_temp);
mimg_c = zeros(dat.ops.Ly, dat.ops.Lx);
mimg_c(dat.ops.yrange, dat.ops.xrange) = squeeze(dat.mimg(:,:,2));

%mimg_c2 = zeros(dat.ops.yrange, dat.ops.xrange);
if isfield(dat.ops, 'scaleXY')
    mimg_c = imresize(mimg_c, 100/dat.ops.scaleXY);
end
if dat.ops.pockelsLineBlank > 0
    croppedX = ceil(0.5*dat.ops.Lx_ori * dat.ops.pockelsLineBlank/100)+1 ...
        :floor(dat.ops.Lx_ori - 0.5*dat.ops.Lx_ori * dat.ops.pockelsLineBlank/100);
    mimg_c2 = zeros(dat.ops.Ly_ori, dat.ops.Lx_ori);
    mimg_c2(:,croppedX) = mimg_c;
    mimg_c = mimg_c2;
end
mimgOriginalSize = mimg_c;