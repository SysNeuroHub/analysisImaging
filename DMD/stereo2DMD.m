%% subject info
subject = 'tmpD';
regDir = '/home/daisuke/Documents/git/analysisImaging/MROIDMD/'; %will be uploaded to market
load(fullfile(regDir, subject, 'Atlas_reg_info.mat'), 'proj_brain','ROI_info',...
    'tform','tform2','mrwarpedtoDMD', 'mrwarped','image2','mrangle');
OIsize = size(image2);
MmPerPixel_oi = 0.0104;

%% image in stereotaxic coordinates (common across subjects)
imgPrefix = 'CCFBL_584x450pix_5x8circle_bk'; %  showAllenCCFBregmaLambda_patches
%imgPrefix = 'natural_400x300pix_2'; %showNatural

imgDir = fullfile('/home/daisuke/tmp/',imgPrefix);
load(fullfile(imgDir, [imgPrefix '_stereo']), 'imageStereo','camImg');
% load(fullfile(imgDir, [imgPrefix '_stereo']), 'imageStereo','imageStereo_CCF','bregma', 'MmPerPixel_img');


%% convert stereo image into image for DMD
[image4DMD, image4OI] = applyStereo2DMD(double(imageStereo)/double(intmax("uint8")), camImg.bregmapix, camImg.MmPerPixel, ...
    mrangle, tform, tform2, OIsize, MmPerPixel_oi, fullfile(regDir, subject));
%12s per image

image4DMD = uint8(round(double(intmax("uint8"))*image4DMD));

%% check if images are within the DMD projection zone
%ideally this check should be done much earlier...
% DMDprojectionZone = getDMDprojectionZone(camImg);
% imagesc(imresize(sum(image4OI,3),.5))
% hold on
% contour(DMDprojectionZone,'r');
% plot(camImg.bregmapix(2), camImg.bregmapix(1), 'yx')
% plot(camImg.lambdapix(2), camImg.lambdapix(1), 'yx');
% %< image4OI: why DMD ptn on the midline does not align bregma and lambda?
% > because tmpD reference image was not aligned to the midline...

save(fullfile(imgDir, [imgPrefix '_' subject]), ...
    'image4DMD','image4OI'); %what else to save?

binary = 1;
mkdir(fullfile(imgDir, subject));
saveEveryImages(image4DMD, fullfile(imgDir, subject), binary, [imgPrefix '_' subject]); %is this really needed?
