function TexImgwarpedtoDMD = applyStereo2DMD(images, bregma, MmPerPixel_img, ...
    usFactor, tform_T2OI, tform_OIDMD, OIsize, MmPerPixel_oi, tmpDir)

MRIdir = 'C:/Documents/git/analysisImaging/MROIDMD';
%OIsize = [1080 1080];  %[y x]
%MmPerPixel_oi = 0.0104;
DMDsize = [500 800]; %[y x]
MmPerPixel_mr = 0.1; %[mm]

if isempty(MmPerPixel_oi)
    MmPerPixel_oi = MmPerPixel_img;
end
if isempty(tmpDIr)
    tmpDir = './tmp/';
end
mkdir(tmpDir);

copyfile(path2T2w, fullfile(tmpDir,'T2w_brain.nii'));
niftiwrite_us(fullfile(tmpDir,'T2w_brain.nii'), usFactor);

%imagesize = [size(images,1) size(images,2)];
nImages = size(images,3);

[tform3, surfDepth, Vus] = registerStereo2CCF(bregma, MmPerPixel_img, usFactor);

TexImgwarpedtoDMD = zeros(DMDsize(1), DMDsize(2), nImages);
for ii = 1:nImages
    disp(['applyStereo2DMD: ' num2str(ii) '/' num2str(nImages)]);
    testimg_reg = imwarp(images(:,:,ii),tform3,'linear', ...
        'OutputView', refDMDbig, 'FillValues',0); %all
    TexVol = paintSurfaceToVolume(surfDepth, testimg_reg, size(Vus));%all
    info = V_info;
    info.Datatype = 'double';
    info.Transform.T(1:3,1:3) = info.Transform.T(1:3,1:3)/usFactor;
    info.PixelDimensions = info.PixelDimensions/usFactor;
    info.ImageSize = size(TexVol);
    niftiwrite(TexVol,'TexVol.nii', info); %all

    cmdStr = [fullfile(MRIdir,'pattern_generation/AtlasTexVol_to_T2.sh') ' ' ...
        'TexVol.nii' ' ' tmpDir];
    system(cmdStr); %output TexVol_T2.nii %all

    %% project back from 3D to 2D
    TexVol_T2 = niftiread('TexVol_T2.nii');

    %% volume mask returns less noisier image than surface mask
    TexImg = getSufraceData2(TexVol_T2, Vus>0);
    %[~, TexImg] = getSurfaceData3(TexVol_T2, 'last', 0);  %slightly faster but noisier

    %% project to OI
    fixedRef  = imref2d(OIsize, MmPerPixel_oi, MmPerPixel_oi);  % example pixel sizes in mm
    movingRef = imref2d(size(TexImg), MmPerPixel_mr/usFactor,  MmPerPixel_mr/usFactor);
    TexImgwarped = imwarp(TexImg,movingRef,tform_T2OI,'linear','OutputView',fixedRef);

    %% project to DMD
    TexImgwarpedtoDMD(:,:,ii) = imwarp(TexImgwarped,tform_OIDMD,'linear','OutputView',imref2d(DMDsize));
end
delete(tmpDir);