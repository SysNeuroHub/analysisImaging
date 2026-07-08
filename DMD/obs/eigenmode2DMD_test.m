%
subject = 'Confucious';
binary = 0;
verbose = 0;
mrangle = [];

addpath(genpath('/home/daisuke/Documents/git/analysisImaging'));
regDir =fullfile('/home/daisuke/Documents/git/analysisImaging/MROIDMD', subject);


Vmode = niftiread('/home/daisuke/Documents/MATLAB/mouseEigenmode/hemisphere_whole/eigenmode_allen100um_hemi=whole_mode=50_frequency=0.407Hz.nii.gz');
Vmode = flip(permute(Vmode,[1 3 2]),3);

% [surfDepth] = getSurfaceData(V, 'last', 0.5);

path2Brain_template =  fullfile(regDir,'Brain_template.nii');

%% from register2Stereo2CCF
usFactor = 1; %5

V = afterniftiread(niftiread(path2Brain_template)>0);
V_info = niftiinfo(path2Brain_template);
V = smooth3(V, 'gaussian', [7 7 7], 3); %important for projection mapping

[xx,yy,zz] = ndgrid(-2:2);
nhood = sqrt(xx.^2 + yy.^2 + zz.^2) <= 2.0;
Verode = imerode(V, nhood);
[surfDepth] = getSurfaceData(Verode, 'last', 0.95);

% Vus = imresize3(V, usFactor, 'linear');
%[surfDepth] = getSurfaceData(Vus, 'last', 0.5);

V_info = V_info;
V_info.Datatype = 'double';
V_info.ImageSize = size(Verode);


%% from Stereo2T2  TexVol_T2 = Stereo2T2(imageStereo, usFactor, Vusinfo, tform3, surfDepth, mrangle)
DMDsize = [500 800]; %[y x]
MmPerPixel_mr = 0.1; %[mm]

cd(regDir);

%% texture mapping from 2D to 3D in Allen CCF (Kim) space
TexVol_CCF = beforeniftiwrite(extractSurfaceValues(surfDepth, Vmode));
% TexVol_CCF = V(surfDepth) .* Vmode;

     if verbose
         subplot(1,3,3); isosurface(TexVol_CCF); axis equal tight;
         f = gcf;
         f.CurrentAxes.YDir = 'Reverse';
         ax = gca; set(gca, 'view', [120 20]);
         title('Projection mapped to upscaled CCF');
     end

    % niftiwrite(flip(TexVol_CCF,3), 'TexVol.nii', Vusinfo);
    Vusinfo4saving = V_info;
    Vusinfo4saving.ImageSize = [V_info.ImageSize(2) V_info.ImageSize(3) V_info.ImageSize(1)];
    niftiwrite(TexVol_CCF, 'TexVol.nii', Vusinfo4saving);

    %% transform from CCF to individual T2
    cmdStr = ['/home/daisuke/Documents/git/analysisImaging/DMD/AtlasTexVol_to_T2.sh' ' ' 'TexVol.nii' ' '  'T2w_brain.nii' ' ' pwd];
    % cmdStr = [fullfile(fileparts(mfilename('fullpath')), 'AtlasTexVol_to_T2.sh') ' ' 'TexVol.nii' ' ' pwd];
    system(cmdStr); 
    %expect: T2w_brain_us.nii, Tem_to_T21Warp.nii.gz, Tem_to_T20GenericAffine.mat
    %output: TexVol_T2.nii 

    %% project back from 3D to 2D
    TexVol_T2 = niftiread('TexVol_T2.nii');
    % if ~isempty(mrangle)
    %     if mrangle(1)~=0
    %         TexVol_T2 = imrotate3(TexVol_T2, mrangle(1), [1 0 0],'linear','crop','fillvalues',nan); %roll
    %     end
    %     if mrangle(2)~=0
    %         TexVol_T2 = imrotate3(TexVol_T2, mrangle(2), [0 1 0],'linear','crop','fillvalues',nan); %pitch
    %     end
    %     if mrangle(3)~=0
    %         TexVol_T2 = imrotate3(TexVol_T2, mrangle(3), [0 0 1],'linear','crop','fillvalues',nan); %yaw
    %     end
    % end
    % TexVol_T2 = afterniftiread(TexVol_T2);

    %% from applyStereo2DMD
TexImg = squeeze(nanmean(TexVol_T2,3));%./surfaceT2,3));
    if verbose
        subplot(1,3,1); imagesc(TexImg);axis equal tight; 
        title('CCF warped to upscaled T2');
    end

    %% project to OI
    if autoTform
        fixedRef  = imref2d(OIsize, MmPerPixel_oi, MmPerPixel_oi);  % example pixel sizes in mm
        movingRef = imref2d(size(TexImg), MmPerPixel_mr/usFactor,  MmPerPixel_mr/usFactor);
        TexImgwarped(:,:,ii) = imwarp(TexImg,movingRef,tform_T2OI_us,'linear','OutputView',fixedRef);
    else
        TexImgwarped(:,:,ii) = imwarp(TexImg,tform_T2OI_us,'linear','OutputView',imref2d(OIsize));
    end

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