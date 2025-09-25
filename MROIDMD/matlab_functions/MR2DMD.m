%need to copy from Chris's result for individual subject
% pills_labels.nii


subjectName = 'tmpB';

MRIdir = '~/Dropbox/GSlab/MRI_alignment/Kim2023_test';
addpath(genpath(MRIdir));
mkdir(fullfile(MRIdir,subjectName));
cd(fullfile(MRIdir,subjectName));
copyfile(fullfile(MRIdir, 'pattern_generation/Allen_annotation_modified.nii'), ...
    fullfile(MRIdir, subjectName,'Allen_annotation_modified.nii'));
copyfile(fullfile(MRIdir, 'pattern_generation/Brain_template.nii'), ...
    fullfile(MRIdir, subjectName,'Brain_template.nii'));
copyfile(fullfile(MRIdir, 'matlab_functions/Allen_pills_mask.nii'), ...
    fullfile(MRIdir, subjectName,'Allen_pills_mask.nii'));

%% subject-specific data
if strcmp(subjectName, 'tmpB')
    nii_ori = '/mnt/dshi0006_market/MRI/analysis/MRA056-selection/20250512-testB/20250512_152043_MRA056_test_20250512B_1_3/8/pdata/1/nifti/MRA056_test_20250512B_8_1_1.nii';
    angle = [0 -(90-79) 0]; %roll pitch yaw
    load('/mnt/dshi0006_market/Subjects/robita/2025-03-29_1/5/dataSummary_amber.mat');
    %J=60;
elseif strcmp(subjectName, 'tmpC')
    nii_ori = '/mnt/dshi0006_market/MRI/analysis/MRA056-selection/20250512-testC/20250512_155921_MRA056_test_20250512C_1_4/6/pdata/1/nifti/MRA056_test_20250512C_6_1_1.nii';
    angle = [0 -(90-85) 0]; %roll pitch yaw
    load('/mnt/dshi0006_market/Subjects/yamatotakeru/2025-04-21_1/1/dataSummary_amber.mat');
    %J=60;
elseif strcmp(subjectName, 'tmpD')
    nii_ori = '/mnt/dshi0006_market/MRI/analysis/MRA056-selection/20250512-testD/20250512_163121_MRA056_test_20250512D_1_5/7/pdata/1/nifti/MRA056_test_20250512D_7_1_1.nii';
    angle = [0 -(90-84) 0]; %roll pitch yaw
    load('/mnt/dshi0006_market/Subjects/mupi/2025-05-05_1/4/dataSummary_amber.mat');
    %J=60;
end

%% prepare Atlas_anno_to_T2.nii & T2w_resample.nii (takes ~5min)
cmdStr = [fullfile(MRIdir,'pattern_generation/Atlas_T2_coreg_DS.sh') ' ' nii_ori ' ' fullfile(MRIdir, subjectName)];
system(cmdStr);

delete(fullfile(MRIdir, subjectName,'Allen_annotation_modified.nii'));
delete(fullfile(MRIdir, subjectName,'Brain_template.nii'));
delete(fullfile(MRIdir, subjectName,'Allen_pills_mask.nii'));

%% load widefield image w reference tubes
image2 = (dataSummary.meanImage);
image2 = image2-min(image2(:));
image2 = image2/max(image2(:));

%% load DMD ref image 
image3 = rgb2gray(imread('star_800x500.png')); %

%% load project DMD ref image captured by widefield camera
image4 = double(imread('star_800x500.tif'));
image4 = image4/max(image4(:));
image4 = image4(:,3:end-2); %hack to align pixel size to image2

%% register MR - OI - DMD
load_mr_bead = niftiread('pills_labels.nii')>0;
load_mr_brain = niftiread('T2w_brain.nii');
load_anno = niftiread('Atlas_anno_to_T2.nii');

DMD_pattern_prep(load_mr_bead, load_mr_brain, load_anno, image2, image3, image4, angle);



%% sanity check
load('Atlas_reg_info.mat');

%% individual brain and CCF warped to individual brain (result of Atlas_T2_coreg_DS3)
ax(1)=subplot(121);imagesc(mrimg_brain); axis equal tight; hold on; contour(proj_anno_cortex~=0,'r'); contour(mrimg_bead~=0,'g'); title('Original MR bead image + warped CCF (r)')
ax(2)=subplot(122);imagesc(proj_anno_cortex); hold on; contour(proj_anno_cortex~=0,'r'); axis equal tight; title('CCF warped to individual brain')
colormap(gray);
screen2png(['Atlas_T2_coreg_' subjectName]  );

%% all warped to input image to DMD
figure('Position',[ 71           1        1850         961]);
ax(1)=subplot(221);imagesc(OIwarpedtoDMD);axis equal tight; grid minor; title('OIwarpedtoDMD');
ax(2)=subplot(222);imshowpair(image3,OIwarpedtoDMD);axis equal tight; grid minor; title('OIwarpedtoDMD(m), image3(g)');
ax(3)=subplot(223);imagesc(mrwarpedtoDMD);axis equal tight; grid minor;
hold on; contour(mrwarpedtoDMD>0,'r'); title('mrwarpedtoDMD');
ax(4)=subplot(224);imagesc(100*borigwarpedtoDMD);axis equal tight; grid minor;
hold on; contour(mapconfwarpedtoDMD,1,'r'); contour(mrwarpedtoDMD>0,'r');
title('borigwarpedtoDMD, mrwarpedtoDMD(r)');
screen2png(['warpedtoDMD_' subjectName]);
