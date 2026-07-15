function TexImgwarped = DMD2OI(TexImgwarpedtoDMD, tform_OIDMD, OIsize)
% TexImgwarped = DMD2OI(TexImgwarpedtoDMD, tform_OIDMD, OIsize)
% rerutns image(s) of DMD pattern in OI space
% 
% INPUT
% TexImgwarpedtoDMD: DMD stimulation image(s)
% tform_OIDMD: "tform2" of Atlas_reg_info.mat
% 
% cf. OI > DMD in applyStereo2DMD
% TexImgwarpedtoDMD(:,:,ii) = imwarp(squeeze(TexImgwarped(:,:,ii)) ,tform_OIDMD,'linear', ...
%        'OutputView',imref2d(DMDsize));

nFrames = size(TexImgwarpedtoDMD,3);
DMDsize(1) = size(TexImgwarpedtoDMD,1);
DMDsize(2) = size(TexImgwarpedtoDMD,2);
TexImgwarped = zeros(OIsize(1),OIsize(2), nFrames, class(TexImgwarpedtoDMD));
RDMD = imref2d(DMDsize);
RTex = imref2d(OIsize);
for ii = 1 : nFrames
    TexImgwarped(:,:,ii) =  imwarp(TexImgwarpedtoDMD(:,:,ii), RDMD, ...
        invert(tform_OIDMD),...
        'linear','OutputView',RTex);
    %'nearest' - NG
    % 'linear'
    % 'cubic'
    % 'bilinear'
    % 'bicubic'
end