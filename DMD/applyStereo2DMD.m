function [TexImgwarpedtoDMD, TexImgwarped] = applyStereo2DMD(images, bregma, MmPerPixel_img, ...
    mrangle, tform_T2OI, tform_OIDMD, OIsize, MmPerPixel_oi, MROIDMDsubjectDir)
% TexImgwarpedtoDMD = applyStereo2DMD(images, bregma, MmPerPixel_img, ...
%     mrangle, tform_T2OI, tform_OIDMD, OIsize, MmPerPixel_oi, MROIDMDsubjectDir)
% transforms images defined in stereotaxic coordinates to DMD space

if nargin < 8
    MmPerPixel_oi = [];
end

verbose = 0;

%OIsize = [1080 1080];  %[y x]
%MmPerPixel_oi = 0.0104;
DMDsize = [500 800]; %[y x]
MmPerPixel_mr = 0.1; %[mm]
usFactor = 5;
path2Brain_template =  '/home/daisuke/Documents/git/analysisImaging/MROIDMD/pattern_generation/Brain_template.nii';

if isempty(MmPerPixel_oi)
    MmPerPixel_oi = MmPerPixel_img;
end

if max(images(:)>1)||min(images(:)<0) 
    error('Input image must be [0 1]');
end
nImages = size(images,3);

[tform3, surfDepth, Vusinfo] = registerStereo2CCF(bregma, MmPerPixel_img, ...
    path2Brain_template, usFactor);

%for AtlasTexVol_to_T2.sh
niftiwrite_us(fullfile(MROIDMDsubjectDir,'T2w_brain.nii'), usFactor);

currentDir = pwd;
cd(MROIDMDsubjectDir);

surfaceT2 = Stereo2T2(ones(size(images(:,:,1))), usFactor, Vusinfo, tform3, surfDepth, mrangle);
%surfaceT2(abs(surfaceT2)<1e-1)=nan;

TexImgwarped = zeros(OIsize(1), OIsize(2), nImages);
TexImgwarpedtoDMD = zeros(DMDsize(1), DMDsize(2), nImages);
for ii = 1:nImages
    disp(['applyStereo2DMD: ' num2str(ii) '/' num2str(nImages)]);

    TexVol_T2 = Stereo2T2(images(:,:,ii), usFactor, Vusinfo, tform3, surfDepth, mrangle);


   TexImg = squeeze(nanmean(TexVol_T2./surfaceT2,3));
    if verbose
        subplot(1,3,1); imagesc(TexImg);axis equal tight; 
        title('CCF warped to upscaled T2');
    end

    %% project to OI
    fixedRef  = imref2d(OIsize, MmPerPixel_oi, MmPerPixel_oi);  % example pixel sizes in mm
    movingRef = imref2d(size(TexImg), MmPerPixel_mr/usFactor,  MmPerPixel_mr/usFactor);
    TexImgwarped(:,:,ii) = imwarp(TexImg,movingRef,tform_T2OI,'linear','OutputView',fixedRef);

    if verbose
        subplot(1,3,2); imagesc(TexImgwarped(:,:,ii) );axis equal tight; 
        title('T2 warped to OI'); 
    end

    %% project to DMD
    TexImgwarpedtoDMD(:,:,ii) = imwarp(TexImgwarped(:,:,ii) ,tform_OIDMD,'linear', ...
        'OutputView',imref2d(DMDsize));
    if verbose
        subplot(1,3,3); imagesc(TexImgwarpedtoDMD);axis equal tight; 
        title('OI warped to DMD');

        fig = gcf;
        caxes(fig, prctile(images(:), [0 100]));
    end
end
delete('TexVol.nii', 'TexVol_T2.nii');
cd(currentDir);