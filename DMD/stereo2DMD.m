%% subject info
subject = 'tmpD';
regDir = '/home/daisuke/Documents/git/analysisImaging/MROIDMD/'; %will be uploaded to market
load(fullfile(regDir, subject, 'Atlas_reg_info.mat'),...
    'tform','tform2','mrwarpedtoDMD', 'mrwarped','image2','mrangle');
OIsize = size(image2);
MmPerPixel_oi = 0.0104;

%% stereo image
% imgDir = '/home/daisuke/Documents/git/analysisImaging/DMD/';
% imwidth = 800;
% imheight = 500;
% MmPerPixel_img = 0.0104 * 9/5; 
% bregma = [5/9*(380-20) imwidth/2+1]+0.5; %[y x]
% 
% figName = ['CCFBL_' num2str(imwidth) 'x' num2str(imheight)];
% imagesStereo = rgb2gray(imread(fullfile(imgDir, [figName '.png'])));    
% imagesStereo = double(imagesStereo/max(imagesStereo(:)));
% 
% imagesStereo = repmat(imagesStereo, [1 1 5]); 

imgDir = '/home/daisuke/tmp/CCFBL_400x300pix_8x7grid/';
% created by showAllenCCFBregmaLambda_patches
load(fullfile(imgDir, 'CCFBL_400x300pix_8x7grid_stereo'), 'imageStereo','bregma', 'MmPerPixel_img');

%% convert stereo image into image for DMD
%imageStereo(:,:,1) = 1;
image4DMD = applyStereo2DMD(imageStereo(:,:,1), bregma, MmPerPixel_img, ...
    mrangle, tform, tform2, OIsize, MmPerPixel_oi, fullfile(regDir, subject));

image4DMD = uint8(round(double(intmax("uint8"))*image4DMD));

mkdir(fullfile(imgDir, subject));
save(fullfile(imgDir, ['CCFBL_400x300pix_8x7grid_' subject]), ...
    'image4DMD'); %what else to save?

saveEveryImages(imageStack, fullfile(imgDir, subject)); %is this really needed?