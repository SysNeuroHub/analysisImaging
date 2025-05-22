addpath(genpath('C:\Users\dshi0006\allenCCF'));

t = Tiff('\\vault-v2.erc.monash.edu.au\MNHS-dshi0006\Subjects\L4GCaMP6s_251\2020-09-21\blue.tif','r');
blue = read(t);
blue_t = blue;
blue_t = fliplr(blue_t');
blue_t = double(blue_t);

StereotaxicInfo_t.yawAngle = -3; %[deg]
StereotaxicInfo_t.rollAngle= 20; %[deg]
StereotaxicInfo_t.pitchAngle = 0; %[deg]
StereotaxicInfo_t.bregmaxpix = 504;%474 % position of bregma in acquired image [pix]
StereotaxicInfo_t.bregmaypix = 315; % position of bregma in acquired image [pix]
MmPerPixel_t = 1e-3*5/0.5*2; %pixel size [mm/pix]

%% option 1: show as stereotaxic coords (Shimaoka Cell Rep / elife)
stereox_i = -6:.01:1;%-5:.01:1;
stereoy_i = -3:.01:6;
blue_t_s = impix2stereo_test(blue_t, StereotaxicInfo_t, MmPerPixel_t, ...
    stereox_i, stereoy_i);
imagesc(stereox_i, stereoy_i, blue_t_s);
axis equal tight
drawTopDownCtx;
title('show as stereotaxic coords');

