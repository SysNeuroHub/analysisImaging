function TexImgwarpedtoDMD = applyStereo2DMD(images, bregma, MmPerPixel_img, ...
    mrangle, tform_T2OI, tform_OIDMD, OIsize, MmPerPixel_oi, MROIDMDsubjectDir)

if nargin < 8
    MmPerPixel_oi = [];
end

verbose = 1;

%OIsize = [1080 1080];  %[y x]
%MmPerPixel_oi = 0.0104;
DMDsize = [500 800]; %[y x]
MmPerPixel_mr = 0.1; %[mm]
usFactor = 5;
path2Brain_template =  '/home/daisuke/Documents/git/analysisImaging/MROIDMD/pattern_generation/Brain_template.nii';

if isempty(MmPerPixel_oi)
    MmPerPixel_oi = MmPerPixel_img;
end

nImages = size(images,3);

[tform3, surfDepth, sizeVus, Vusinfo] = registerStereo2CCF(bregma, MmPerPixel_img, ...
    path2Brain_template, usFactor);

%for AtlasTexVol_to_T2.sh
niftiwrite_us(fullfile(MROIDMDsubjectDir,'T2w_brain.nii'), usFactor);

currentDir = pwd;
cd(MROIDMDsubjectDir);

TexImgwarpedtoDMD = zeros(DMDsize(1), DMDsize(2), nImages);
for ii = 1:nImages
    disp(['applyStereo2DMD: ' num2str(ii) '/' num2str(nImages)]);

    if verbose
        figure;
        subplot(2,3,1); imagesc(images(:,:,ii)); hold on; plot(bregma(2), bregma(1),'ro'); axis equal tight;
        title('image in stereotaxic coords');
    end

    %% warp from stereotaxic to CCF (2D)
    refDMD = imref2d(round([sizeVus(3), sizeVus(1)]/usFactor));
    refDMDbig = imref2d([sizeVus(3), sizeVus(1)]);
    refDMDbig.XWorldLimits = refDMD.XWorldLimits;
    refDMDbig.YWorldLimits = refDMD.YWorldLimits;

     img_reg = imwarp(images(:,:,ii),tform3,'linear', ...
        'OutputView', refDMDbig, 'FillValues',0); %all

     if verbose
         subplot(2,3,2); imagesc(img_reg); axis equal tight; 
         title('stereo image registered to upscaled CCF');
     end

    %% texture mapping from 2D to 3D in Allen CCF (Kim) space
     TexVol_CCF = paintSurfaceToVolume(surfDepth, img_reg, sizeVus); %CORRECT
     %TexVol_CCF = imdilate(TexVol_CCF, strel('sphere', 1)); %does not reduce noise in TexImg

     if verbose
         subplot(2,3,3); isosurface(TexVol_CCF); axis equal tight;
         f = gcf;
         f.CurrentAxes.ZDir = 'Reverse';
         ax = gca; set(gca, 'view', [120 20]);
         title('Projection mapped to upscaled CCF');
     end

    niftiwrite(flip(TexVol_CCF,3), 'TexVol.nii', Vusinfo);

    cmdStr = [fullfile(fileparts(mfilename('fullpath')), 'AtlasTexVol_to_T2.sh') ' ' 'TexVol.nii' ' ' pwd];
    system(cmdStr); 
    %expect: T2w_brain_us.nii, Tem_to_T21Warp.nii.gz, Tem_to_T20GenericAffine.mat
    %output: TexVol_T2.nii 

    %% project back from 3D to 2D
    TexVol_T2 = niftiread('TexVol_T2.nii'); %WRONG

    if mrangle(1)~=0
        TexVol_T2 = imrotate3(TexVol_T2, mrangle(1), [1 0 0],'linear','crop'); %roll
    end
    if mrangle(2)~=0
        TexVol_T2 = imrotate3(TexVol_T2, mrangle(2), [0 1 0],'linear','crop'); %pitch
    end
    if mrangle(3)~=0
        TexVol_T2 = imrotate3(TexVol_T2, mrangle(3), [0 0 1],'linear','crop'); %yaw
    end

   %[~, TexImg] = getSurfaceData3(TexVol_T2, 'last', 0.5); %noisy
  TexImg = flipud(squeeze(max(TexVol_T2,[],2))');
   % TexVol_T2(TexVol_T2==0) = NaN;
   % TexImg = flipud(squeeze(nanmean(TexVol_T2,2))');
    if verbose
        subplot(2,3,4); imagesc(TexImg);axis equal tight; 
        title('CCF warped to upscaled T2');
    end

    %% project to OI
    fixedRef  = imref2d(OIsize, MmPerPixel_oi, MmPerPixel_oi);  % example pixel sizes in mm
    movingRef = imref2d(size(TexImg), MmPerPixel_mr/usFactor,  MmPerPixel_mr/usFactor);
    TexImgwarped = imwarp(TexImg,movingRef,tform_T2OI,'linear','OutputView',fixedRef);

    if verbose
        subplot(2,3,5); imagesc(TexImgwarped);axis equal tight; 
        title('T2 warped to OI'); 
    end

    %% project to DMD
    TexImgwarpedtoDMD(:,:,ii) = imwarp(TexImgwarped,tform_OIDMD,'linear', ...
        'OutputView',imref2d(DMDsize));
    if verbose
        subplot(2,3,6); imagesc(TexImgwarpedtoDMD);axis equal tight; 
        title('OI warped to DMD');
    end
end
delete('TexVol.nii', 'TexVol_T2.nii');
cd(currentDir);