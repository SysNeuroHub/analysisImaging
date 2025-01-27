function [image_stereo, xaxis_stereo, yaxis_stereo] = impix2stereo_test2(image, StereotaxicInfo, MmPerPixel)
% image_stereo_i = impix2stereo_test(image, StereotaxicInfo, MmPerPix, stereox_i, stereoy_i)
% returns image aligned 
% inputs:
% StereotaxicInfo.yawAngle, 
% StereotaxicInfo.rollAngle,
% StereotaxicInfo.pitchAngle
% StereotaxicInfo.bregmaxpix, StereotaxicInfo.bregmaypix
% 
% 2016-1-13 DS created
% TODO: rotate pitch roll and yaw at one time

%% rotate brain along yaw
roll = StereotaxicInfo.rollAngle;%23;%deg
pitch = StereotaxicInfo.pitchAngle;%5;%deg
yaw = StereotaxicInfo.yawAngle;
%maskInfo.StereotaxicInfo.yawAngle = -40;%test

image_yaw = imrotate(image, yaw, 'nearest', 'crop');

%yawAngle = StereotaxicInfo.yawAngle*pi/180;

%% rotate bregma position along yaw
%this is only applicable if imrotate is done with crop option
sz = size(image_yaw) / 2;
orig_x = StereotaxicInfo.bregmaxpix - sz(2);
orig_y = StereotaxicInfo.bregmaypix - sz(1);
rot_mat=[cosd(-yaw), sind(-yaw); -sind(-yaw) ,cosd(-yaw)];
old_orig = [orig_x orig_y];
new_orig =old_orig * rot_mat;
bregmaxpix_yaw = new_orig(1) + sz(2);
bregmaypix_yaw = new_orig(2) + sz(1);


xaxis_yaw = ((1:size(image_yaw,2)) - bregmaxpix_yaw) * MmPerPixel; %mm from bregma
yaxis_yaw = ((1:size(image_yaw,1)) - bregmaypix_yaw) * MmPerPixel; %mm from bregma



%% rotate in pitch and roll
tform = projective2d([cosd(roll) 0 0;0 cosd(pitch) 0; 0 0 1]);%this may be inaccurate
coordinate_yaw = imref2d(size(image_yaw), ...
    [xaxis_yaw(1) xaxis_yaw(end)], [yaxis_yaw(1) yaxis_yaw(end)]);
[image_stereo, coordinate_stereo] = imwarp(image_yaw, coordinate_yaw, tform);
[bregmaxpix_stereo, bregmaypix_stereo] = coordinate_stereo.worldToIntrinsic(0,0);
xaxis_stereo = ((1:size(image_stereo,2)) - bregmaxpix_stereo) * MmPerPixel; %mm from bregma. pixel size does not change after rotation
yaxis_stereo = ((1:size(image_stereo,1)) - bregmaypix_stereo) * MmPerPixel; %mm from bregma. pixel size does not change after rotation


% %% stereotaxic coordinates of the brain image
% stereox = -xaxis_stereo;
% stereoy = -yaxis_stereo;
% [STEREOX, STEREOY] = meshgrid(stereox, stereoy);
% 
% 
% %% interpolate in 2D
% [STEREOX_i, STEREOY_i] = meshgrid(stereox_i, stereoy_i);
% 
% image_stereo_i = interp2(STEREOX, STEREOY, image_stereo, STEREOX_i, STEREOY_i, ...
%     interpMethod, 0);%changed from nearest to spline 25/1/17

% %sanity check
% subplot(221);
% imagesc(image);hold on
% plot(StereotaxicInfo.bregmaxpix, StereotaxicInfo.bregmaypix,'ro');
% xlim([1 size(image,2)+50]);
% axis equal tight;
% title('original')
% 
% subplot(222);
% imagesc(xaxis_yaw, yaxis_yaw,image_yaw);
% hold on
% plot(0,0,'ro')
% axis equal tight
% title('rotate yaw')
% 
% subplot(223);
% imagesc(xaxis_stereo, yaxis_stereo,image_stereo);
% axis equal tight
% hold on
% plot(0,0,'ro')
% title('rotate pitch roll')
% 
% subplot(224);
% imagesc(stereox_i, stereoy_i, image_stereo_i);
% axis equal tight xy


%%NG
% ResizeFac = 1;
% [pixx, pixy] = stereoTaxicCoords(Exps, STEREOX, STEREOY, FileString, ResizeFac);


