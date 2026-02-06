%% subject info
subject = 'tmpD';
regDir = '/home/daisuke/Documents/git/analysisImaging/MROIDMD/'; %will be uploaded to market
load(fullfile(regDir, subject, 'Atlas_reg_info.mat'), 'proj_brain','ROI_info',...
    'tform','tform2','mrwarpedtoDMD', 'mrwarped','image2','mrangle');
OIsize = size(image2);
MmPerPixel_oi = 0.0104;

%% stereo image
imgPrefix = 'CCFBL_584x450pix_9x8circle'; %  showAllenCCFBregmaLambda_patches
%imgPrefix = 'natural_400x300pix_2'; %showNatural

imgDir = fullfile('/home/daisuke/tmp/',imgPrefix);
load(fullfile(imgDir, [imgPrefix '_stereo']), 'imageStereo','camImg');
% load(fullfile(imgDir, [imgPrefix '_stereo']), 'imageStereo','imageStereo_CCF','bregma', 'MmPerPixel_img');


%% convert stereo image into image for DMD
[image4DMD, image4OI] = applyStereo2DMD(double(imageStereo)/double(intmax("uint8")), camImg.bregmapix, camImg.MmPerPixel, ...
    mrangle, tform, tform2, OIsize, MmPerPixel_oi, fullfile(regDir, subject));
%12s per image

image4DMD = uint8(round(double(intmax("uint8"))*image4DMD));

mkdir(fullfile(imgDir, subject));
save(fullfile(imgDir, [imgPrefix '_' subject]), ...
    'image4DMD','image4OI'); %what else to save?

binary = 1;
saveEveryImages(image4DMD, fullfile(imgDir, subject), binary); %is this really needed?

%% convert stereo CCF image as a reference
[image4DMD_CCF, image4OI_CCF] = applyStereo2DMD(imageStereo_CCF, bregma, MmPerPixel_img, ...
    mrangle, tform, tform2, OIsize, MmPerPixel_oi, fullfile(regDir, subject));
image4DMD_CCF = uint8(round(double(intmax("uint8"))*image4DMD_CCF));
save(fullfile(imgDir, [imgPrefix '_' subject]), ...
    'image4DMD_CCF','image4OI_CCF', '-append');

