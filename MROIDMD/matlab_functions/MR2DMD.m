%need to copy from Chris's result for individual subject
% pills_labels.nii



subjectName = 'tmpD';
if ispc
    MRdataServer = 'M:/';
else
    MRdataServer = '/mnt/dshi0006_market';
    MRIdir = '/home/daisuke/Documents/git/analysisImaging/MROIDMD';
    refDir = '/home/daisuke/Documents/git/analysisImaging/DMD/references';
end
refdate = '20260214';

%% camera image info
load(fullfile(refDir, ['camImg_' refdate]),'camImg');

addpath(genpath(MRIdir));
mkdir(fullfile(MRIdir,subjectName));
cd(fullfile(MRIdir,subjectName));
copyfile(fullfile(MRIdir, 'pattern_generation/Allen_annotation_modified.nii'), ...
    fullfile(MRIdir, subjectName,'Allen_annotation_modified.nii'));
copyfile(fullfile(MRIdir, 'pattern_generation/Brain_template.nii'), ...
    fullfile(MRIdir, subjectName,'Brain_template.nii'));
copyfile(fullfile(MRIdir, 'matlab_functions/Allen_pills_mask.nii'), ...
    fullfile(MRIdir, subjectName,'Allen_pills_mask.nii'));

%% path setting for .sh scripts
%setenv('PATH', [getenv('PATH')]);
setenv('PATH',  '/home/daisuke/anaconda3/bin:/home/daisuke/anaconda3/condabin:/usr/local/freesurfer/8.0.0/bin:/usr/local/freesurfer/8.0.0/fsfast/bin:/usr/local/freesurfer/8.0.0/tktools:/home/daisuke/fsl/share/fsl/bin:/usr/local/freesurfer/8.0.0/mni/bin:/home/daisuke/fsl/share/fsl/bin:/home/daisuke/fsl/share/fsl/bin:/home/daisuke/.local/bin:/home/daisuke/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:/snap/bin:/snap/bin:/home/daisuke/fsl/bin :~/bin :/home/daisuke/abin')
%check
system('which 3dcopy');

%% subject-specific data
if strcmp(subjectName, 'tmpB')
    nii_ori = fullfile(MRdataServer, 'MRI/analysis/MRA056-selection/20250512-testB/20250512_152043_MRA056_test_20250512B_1_3/8/pdata/1/nifti/MRA056_test_20250512B_8_1_1.nii');
    mrangle = [0 -(90-79) 0]; %roll pitch yaw
    %load('/mnt/dshi0006_market/Subjects/robita/2025-03-29_1/5/dataSummary_amber.mat');
    image2 = double(imread('/mnt/dshi0006_vault/Subjects/robita/2025-10-15/redamber_1168x900.tif'));
    image2 = image2-min(image2(:));
    image2 = image2/max(image2(:));
    %J=60;
elseif strcmp(subjectName, 'tmpC')
    nii_ori = fullfile(MRdataServer,'MRI/analysis/MRA056-selection/20250512-testC/20250512_155921_MRA056_test_20250512C_1_4/6/pdata/1/nifti/MRA056_test_20250512C_6_1_1.nii');
    mrangle = [0 -(90-85) 0]; %roll pitch yaw
    %load('/mnt/dshi0006_market/Subjects/yamatotakeru/2025-04-21_1/1/dataSummary_amber.mat');
    image2 = double(imread('/mnt/dshi0006_vault/Subjects/yamatotakeru/2025-10-15/redamber_1168x900.tif'));
    image2 = image2-min(image2(:));
    image2 = image2/max(image2(:));
    %J=60;
elseif strcmp(subjectName, 'tmpD')
    nii_ori = fullfile(MRdataServer,'MRI/analysis/MRA056-selection/20250512-testD/20250512_163121_MRA056_test_20250512D_1_5/7/pdata/1/nifti/MRA056_test_20250512D_7_1_1.nii');
    mrangle = [0 -(90-84) 0]; %roll pitch yaw
    
    %% load widefield image w reference tubes
    % image2 = double(imread('/mnt/dshi0006_vault/Subjects/mupi/2025-10-02/a_24194.TIF'));
    image2 = double(imread('/mnt/dshi0006_vault/Subjects/mupi/2025-10-15/amber_1168x900.TIF'));
    image2 = image2-min(image2(:));
    image2 = image2/max(image2(:));
    %load('/mnt/dshi0006_market/Subjects/mupi/2025-05-05_1/4/dataSummary_amber.mat');
    %J=60;
elseif strcmp(subjectName,'Nero')
    nii_ori = fullfile(MRdataServer, 'MRI/record/20251001-Nero/20251001_141337_MRA056_Nero_20251001_1_7/11/pdata/1/nifti/MRA056_Nero_20251001_11_1_1.nii');
    mrangle = [];
    %dataSummary.meanImage = imread('/mnt/dshi0006_vault/Subjects/Nero/2025-10-01_1/amber.tif');
    dataSummary.meanImage = imread('/mnt/dshi0006_vault/Subjects/Nero/2025-10-02/1000x1600.TIF');
elseif strcmp(subjectName, 'WT78')
    %nii_ori = fullfile(dataServer, '/MRI/record/20260116_151643_MRA056_WT78_20260116_1_8/21/pdata/1/nifti/MRA056_WT78_20260116_21_1_1.nii');
    nii_ori = fullfile(MRdataServer, '/MRI/record/20260116_151643_MRA056_WT78_20260116_1_8/21/pdata/1/nifti/MRA056_WT78_20260116_21_1_1.nii');
    mrangle = [];
    image2 = double(imread('/mnt/dshi0006_vault/Subjects/WT78/amber_1184x900.TIF'));
    image2 = image2-min(image2(:));
    image2 = image2/max(image2(:));
elseif strcmp(subjectName, 'WT79')
    %nii_ori = fullfile(dataServer, '/MRI/record/20260116_151643_MRA056_WT78_20260116_1_8/21/pdata/1/nifti/MRA056_WT78_20260116_21_1_1.nii');
    nii_ori = fullfile(MRdataServer, '/MRI/record/20260116_165900_MRA056_WT79_20260116_1_9/5/pdata/1/nifti/MRA056_WT79_20260116_5_1_1.nii');
    mrangle = [];
    image2 = double(imread('/mnt/dshi0006_vault/Subjects/WT79/amber_1184x900.TIF'));
    image2 = image2-min(image2(:));
    image2 = image2/max(image2(:));
end

%% load project DMD ref image captured by widefield camera
% image4 = double(imread('star_1080x1080.tif')); %dimension must be same as image2
% image4 = double(imread('~/Documents/git/analysisImaging/MROIDMD/matlab_functions/star_1168x900.tif')); %dimension must be same as image2
%image4 = double(imread('~/Documents/git/analysisImaging/MROIDMD/matlab_functions/star_800x500_1168x900.tif')); %22/10/25
% image4 = double(imread('~/Documents/git/analysisImaging/MROIDMD/matlab_functions/star_800x500_1168x900_20251027.tif')); %27/10/25
%image4 = double(imread('~/Documents/git/analysisImaging/MROIDMD/matlab_functions/star_800x500_1168x900_20251028.tif')); %28/10/25
image4_tmp = double(imread(fullfile(refDir, ['star_' refdate '.tif'])));
image4_tmp = resizeCropFullImg(image4_tmp, camImg);
image4_tmp = image4_tmp/max(image4_tmp(:));

%% prepare Atlas_anno_to_T2.nii & T2w_resample.nii (takes ~5min)
cmdStr = [fullfile(MRIdir,'pattern_generation/Atlas_T2_coreg_DS.sh') ' ' nii_ori ' ' fullfile(MRIdir, subjectName)];


system(cmdStr); %calls FindPillsExp_Allen.py inside

delete(fullfile(MRIdir, subjectName,'Allen_annotation_modified.nii'));
delete(fullfile(MRIdir, subjectName,'Brain_template.nii'));
delete(fullfile(MRIdir, subjectName,'Allen_pills_mask.nii'));



%% load DMD ref image 
image3 = rgb2gray(imread(fullfile(MRIdir, 'matlab_functions','star_800x500.png')); %fixed

 
%% register MR - OI - DMD
load_mr_bead = niftiread('pills_labels.nii')>0;
load_mr_brain = niftiread('T2w_brain.nii');
load_anno = niftiread('Atlas_anno_to_T2.nii');

%Atlas_reg_info = 
DMD_pattern_prep(load_mr_bead, load_mr_brain, load_anno, image2, image3, image4, mrangle);
%save Atlas_reg_info.mat


%% sanity check
load('Atlas_reg_info.mat');

%% individual brain and CCF warped to individual brain (result of Atlas_T2_coreg_DS3)
ax(1)=subplot(121);imagesc(mrimg_brain); axis equal tight; hold on; contour(proj_anno_cortex~=0,'r'); contour(mrimg_bead~=0,'g'); title('Original MR bead image + warped CCF (r)')
ax(2)=subplot(122);imagesc(proj_anno_cortex); hold on; contour(proj_anno_cortex~=0,'r'); axis equal tight; title('CCF warped to individual brain')
colormap(gray);
screen2png(['Atlas_T2_coreg_' subjectName]  );

%% all warped to 2D OI space
mr_brain = niftiread('T2w_brain.nii');
mr_brain = imrotate3(mr_brain, mrangle(1), [1 0 0],'linear','crop'); %roll
mr_brain = imrotate3(mr_brain, mrangle(2), [0 1 0],'linear','crop'); %pitch
mr_brain = imrotate3(mr_brain, mrangle(3), [0 0 1],'linear','crop'); %yaw
[mrimg_surf] = getSurfaceData(afterniftiread(mr_brain), 'last');

    fixedRef  = imref2d(size(image2), 0.0104, 0.0104);  % example pixel sizes in mm
movingRef = imref2d(size(mrimg_brain), 0.1, 0.1);
beadwarped = imwarp(mrimg_bead,movingRef,tform,'cubic','OutputView',fixedRef);
surfwarped = imwarp(mrimg_surf,movingRef,tform,'cubic','OutputView',fixedRef);
figure('position', [675         453        1240         413]);
subplot(121);
imagesc(image2);axis equal tight; hold on;
clim(prctile(image2(:), [1 97]));
contour(surfwarped>0, 'y'); contour(beadwarped>.5,'r');
title('widefield image + brain in MR(y) + beads in MR(r)')
subplot(122);
imagesc(surfwarped+beadwarped);
axis equal tight; hold on;
contour(surfwarped>0, 'y'); contour(beadwarped>.5,'r');
title('brain in MR(y) + beads in MR(r)');
colormap(gray);
screen2png(['warpedtoOI_' subjectName]);


%% all warped to input image to DMD
figure('Position',[ 71           1        1850         961]);
ax(1)=subplot(221);imagesc(OIwarpedtoDMD);axis equal tight; grid minor; title('OIwarpedtoDMD');
ax(2)=subplot(222);imshowpair(image3,OIwarpedtoDMD);axis equal tight; grid minor; title('OIwarpedtoDMD(m), image3(g)');
ax(3)=subplot(223);imagesc(mrwarpedtoDMD);axis equal tight; grid minor;
hold on; contour(mrwarpedtoDMD>0,'r'); title('mrwarpedtoDMD');
ax(4)=subplot(224);imagesc(100*borigwarpedtoDMD);axis equal tight; grid minor;
hold on; contour(mapconfwarpedtoDMD,1,'g'); contour(mrwarpedtoDMD>0,'r');
title('borigwarpedtoDMD, mrwarpedtoDMD(r)');
screen2png(['warpedtoDMD_' subjectName]);


% %% check stimuli after running stereo2DMD
% imgPrefix = 'CCFBL_400x300pix_8x7grid';
% imgDir = fullfile('/home/daisuke/tmp/',imgPrefix);
% load(fullfile(imgDir, [imgPrefix '_' subjectName]), ...
%     'image4DMD_CCF','image4OI_CCF');
% figure('position', [675         453        1240         413]);
% subplot(121);
% imagesc(image2);axis equal tight; hold on;
% clim(prctile(image2(:), [1 97]));
% contour(image4OI_CCF>.5, 'm');
% colormap(gray);
% screen2png(['CCFwarpedtoOI_' subjectName]);
