function TexVol_T2 = Stereo2T2(imageStereo, usFactor, Vusinfo, tform3, surfDepth, mrangle)

verbose = 0;
sizeVus = Vusinfo.ImageSize;

   if verbose
        fig = figure;
        subplot(1,3,1); imagesc(imageStereo); hold on; axis equal tight;
        title('image in stereotaxic coords');
    end

    %% warp from stereotaxic to CCF (2D)
    refDMD = imref2d(round([sizeVus(1), sizeVus(2)]/usFactor));
    refDMDbig = imref2d([sizeVus(1), sizeVus(2)]);
    refDMDbig.XWorldLimits = refDMD.XWorldLimits;
    refDMDbig.YWorldLimits = refDMD.YWorldLimits;

     img_reg = imwarp(imageStereo,tform3,'linear', ...
        'OutputView', refDMDbig, 'FillValues',0);

     if verbose
         subplot(1,3,2); imagesc(img_reg); axis equal tight; 
         title('stereo image registered to upscaled CCF');
     end

    %% texture mapping from 2D to 3D in Allen CCF (Kim) space
     TexVol_CCF = beforeniftiwrite(paintSurfaceToVolume(surfDepth, img_reg, sizeVus)); 
     
     if verbose
         subplot(1,3,3); isosurface(TexVol_CCF); axis equal tight;
         f = gcf;
         f.CurrentAxes.YDir = 'Reverse';
         ax = gca; set(gca, 'view', [120 20]);
         title('Projection mapped to upscaled CCF');
     end

    % niftiwrite(flip(TexVol_CCF,3), 'TexVol.nii', Vusinfo);
    Vusinfo4saving = Vusinfo;
    Vusinfo4saving.ImageSize = [Vusinfo.ImageSize(2) Vusinfo.ImageSize(3) Vusinfo.ImageSize(1)];
    niftiwrite(TexVol_CCF, 'TexVol.nii', Vusinfo4saving);

    %% transform from CCF to individual T2
    cmdStr = [fullfile(fileparts(mfilename('fullpath')), 'AtlasTexVol_to_T2.sh') ' ' 'TexVol.nii' ' ' pwd];
    system(cmdStr); 
    %expect: T2w_brain_us.nii, Tem_to_T21Warp.nii.gz, Tem_to_T20GenericAffine.mat
    %output: TexVol_T2.nii 

    %% project back from 3D to 2D
    TexVol_T2 = niftiread('TexVol_T2.nii');
    if mrangle(1)~=0
        TexVol_T2 = imrotate3(TexVol_T2, mrangle(1), [1 0 0],'linear','crop','fillvalues',nan); %roll
    end
    if mrangle(2)~=0
        TexVol_T2 = imrotate3(TexVol_T2, mrangle(2), [0 1 0],'linear','crop','fillvalues',nan); %pitch
    end
    if mrangle(3)~=0
        TexVol_T2 = imrotate3(TexVol_T2, mrangle(3), [0 0 1],'linear','crop','fillvalues',nan); %yaw
    end
    
    TexVol_T2 = afterniftiread(TexVol_T2);