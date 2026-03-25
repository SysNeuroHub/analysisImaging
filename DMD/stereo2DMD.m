%
resultServer = '/mnt/dshi0006_market';


%% subject info
subject = 'Confucious';
binary = 1;

regDir =fullfile('/home/daisuke/Documents/git/analysisImaging/MROIDMD', subject);
% regDir = fullfile('/mnt/dshi0006_market/Subjects',subject,'MR2DMDresult');%better to be somewhere as applyStereo2DMD needs to cd to this directory
load(fullfile(regDir, 'Atlas_reg_info.mat'), 'proj_brain','ROI_info',...
    'tform','tform2','mrwarpedtoDMD', 'mrwarped','image2','mrangle','autoTform','image2th');
OIsize = size(image2);

%% image in stereotaxic coordinates (common across subjects)
% imgPrefix = 'CCFBL_584x450pix_1x3circle_l'; %  showAllenCCFBregmaLambda_patches
imgPrefix = 'CCFBL_584x450pix_1x1circle_r'; %  showAllenCCFBregmaLambda_patches
%imgPrefix = 'natural_400x300pix_2'; %showNatural

% imgDir = fullfile('/mnt/dshi0006_market/DMD images/',imgPrefix); %unstable in reading nifti ... better to use a local directory
imgDir = fullfile('~/tmp',imgPrefix);
load(fullfile(imgDir, [imgPrefix '_stereo']), 'imageStereo','camImg');
MmPerPixel_oi = camImg.MmPerPixel;


%% convert stereo image into image for DMD
[image4DMD, image4OI] = applyStereo2DMD(double(imageStereo./max(imageStereo(:))), camImg.bregmapix, camImg.MmPerPixel, ...
    mrangle, tform, tform2, OIsize, MmPerPixel_oi, regDir, autoTform);
%12s per image


%% check if images are within the DMD projection zone ... ideally this check should be done much earlier...
all_w_CCF = imread(fullfile(imgDir, [imgPrefix '_stereo_all_wCCF.png']));
all_w_CCF(1,:) = 0; all_w_CCF(end,:) = 0; %HACK

[~,image4OI_all_wCCF] = applyStereo2DMD(double(all_w_CCF./max(all_w_CCF(:))), camImg.bregmapix, camImg.MmPerPixel, ...
    mrangle, tform, tform2, OIsize, MmPerPixel_oi, regDir, autoTform);

f = figure;
DMDprojectionZone = getDMDprojectionZone(camImg);
imagesc(normalize_prctile(image2,[image2th,100])); colormap(gray); axis equal tight; hold on
contour(sum(image4OI_all_wCCF, 3),[.5 .5],'w')
contourf(sum(image4OI,3),[.5 5],'c', 'linewidth',1)
contour(DMDprojectionZone,'r');
title('predicted DMD images (c) and projection zone (r) on OI');
plot(camImg.bregmapix(2), camImg.bregmapix(1), 'yx')
% plot(camImg.lambdapix(2), camImg.lambdapix(1), 'yx'); %DO NOT USE
exportPng4DMD(fullfile(imgDir, [imgPrefix '_' subject '_all_wCCF']), f, 0);

f = figure;
DMDprojectionZone = getDMDprojectionZone(camImg);
%imagesc(log(normalize_prctile(image2,[20 image2th]))); colormap(gray); axis equal tight; hold on
imagesc(normalize_prctile(image2,[20, image2th])); colormap(gray); axis equal tight; hold on
 contour(sum(image4OI_all_wCCF, 3),[.5 .5],'w')
 contourf(sum(image4OI,3),[.5 5],'c', 'linewidth',1)
contour(DMDprojectionZone,'r');
title('predicted DMD images (c) and projection zone (r) on OI');
plot(camImg.bregmapix(2), camImg.bregmapix(1), 'yx')
% plot(camImg.lambdapix(2), camImg.lambdapix(1), 'yx'); %DO NOT USE
exportPng4DMD(fullfile(imgDir, [imgPrefix '_' subject '_all_wCCF_brain']), f, 0);

%% save binarized png images for DMD
save(fullfile(imgDir, [imgPrefix '_' subject]), ...
    'image4DMD','image4OI','image4OI_all_wCCF', 'camImg'); %what else to save?

mkdir(fullfile(imgDir, subject));
saveEveryImages(image4DMD, fullfile(imgDir, subject), binary, [imgPrefix '_' subject]); 


%% upload results to Market
copyfile(imgDir, fullfile(resultServer, 'DMD images',imgPrefix));

% mkdir(fullfile(resultServer, 'DMD images',imgPrefix,subject));
% copyfile(fullfile(fullfile(imgDir, subject)), fullfile(resultServer, 'DMD images',imgPrefix,subject));