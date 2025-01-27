
load('sample left hem roll20deg', 'brainImage');

StereotaxicInfo.yawAngle = -5; %[deg]
StereotaxicInfo.rollAngle = 20; %[deg]
StereotaxicInfo.pitchAngle = 0; %[deg]
StereotaxicInfo.bregmaxpix = 485; % position of bregma in acquired image [pix]
StereotaxicInfo.bregmaypix = 315; % position of bregma in acquired image [pix]
MmPerPixel_t = 1e-3*5/0.5*2; %pixel size [mm/pix]

stereox_i = -6:.01:1; % range of stereotaxic coordinate in L-M [mm]
stereoy_i = -5:.01:6; % range of stereotaxic coordinate in A-P [mm]

% convert image to stereotaxic coordinate
brainImage_stereo = impix2stereo_test(brainImage, StereotaxicInfo, MmPerPixel_t, ...
    stereox_i, stereoy_i);
imagesc(stereox_i, stereoy_i, brainImage_stereo);
axis equal tight

%% superimpose CCF
drawTopDownCtx;
plot(0,0, 'ro');%bregma
%axis off