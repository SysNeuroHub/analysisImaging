
stereox_i = -8:.01:1;%-5:.01:1;
stereoy_i = -3:.01:8;%-3:.01:6;


%tilted 20deg
t = Tiff('\\vault-v2.erc.monash.edu.au\MNHS-dshi0006\Subjects\L4GCaMP6s_250\2020-09-21\blue nocone.tif','r');
blue = read(t);
blue_t = blue;
blue_t = fliplr(blue_t');
blue_t = double(blue_t);
subplot(231)
imagesc(blue_t); axis equal tight;
title('original image, tilted');

StereotaxicInfo_t.yawAngle = -5; 
StereotaxicInfo_t.rollAngle= 20;
StereotaxicInfo_t.pitchAngle = 0;
StereotaxicInfo_t.bregmaxpix = 485;
StereotaxicInfo_t.bregmaypix = 315;
MmPerPixel_t = 1e-3*5/0.5*2;

blue_t_s = impix2stereo_test(blue_t, StereotaxicInfo_t, MmPerPixel_t, ...
    stereox_i, stereoy_i);
subplot(234)
imagesc(stereox_i, stereoy_i, blue_t_s);
axis equal tight
grid on
title('converted image in stereo coords');


%% UNtilted
t = Tiff('\\vault-v2.erc.monash.edu.au\MNHS-dshi0006\Subjects\L4GCaMP6s_250\2020-10-9\blue nocone.tif','r');
blue = read(t);
blue = fliplr(blue');
blue = double(blue);
subplot(232)
imagesc(blue); axis equal tight;
title('original image, UNtilted');

StereotaxicInfo.yawAngle = -4; 
StereotaxicInfo.rollAngle= 0;
StereotaxicInfo.pitchAngle = 0;
StereotaxicInfo.bregmaxpix = 813;
StereotaxicInfo.bregmaypix = 825;
MmPerPixel = 1e-3*5/0.5*1;

% show as original pixel resolution
% [blue_s, xaxis_stereo, yaxis_stereo] = impix2stereo_test2(blue, StereotaxicInfo, MmPerPixel);
% imagesc(xaxis_stereo, yaxis_stereo, blue_s);

%show as in stereotaxic
blue_s = impix2stereo_test(blue, StereotaxicInfo, MmPerPixel, ...
    stereox_i, stereoy_i);
subplot(235);
imagesc(stereox_i, stereoy_i, blue_s);
axis equal tight
grid on

subplot(236);
imshowpair(blue_s, blue_t_s);
title('superposition of two images in stereo coords');

 screen2png('alignTwoImages');