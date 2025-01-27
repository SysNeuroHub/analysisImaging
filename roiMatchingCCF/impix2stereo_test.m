function image_stereo_i = impix2stereo_test(image, StereotaxicInfo, MmPerPixel, ...
    stereox_i, stereoy_i, interpMethod)
% image_stereo_i = impix2stereo_test(image, StereotaxicInfo, MmPerPix, stereox_i, stereoy_i)
% returns image in stereotaxic coordinate specified as (stereox_i,stereoy_i)
% inputs:
% image: image of brain surface (A-P, L-M)
% StereotaxicInfo.yawAngle: [deg]
% StereotaxicInfo.rollAngle: [deg]
% StereotaxicInfo.pitchAngle: [deg]
% StereotaxicInfo.bregmaxpix: location(x) of bregma in image [pixel]
% StereotaxicInfo.bregmaypix: location(y) of bregma in image [pixel]
% MmPerPixel: pixel size of image [mm]
%
% 2016-1-13 DS created

if nargin < 6
    interpMethod = 'spline';
end

%% rotate brain along yaw
roll = StereotaxicInfo.rollAngle;
pitch = StereotaxicInfo.pitchAngle;
yaw = StereotaxicInfo.yawAngle;

image_yaw = imrotate(image, yaw, 'nearest', 'crop');


%% rotate bregma position along yaw
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
xaxis_stereo = ((1:size(image_stereo,2)) - bregmaxpix_stereo) * MmPerPixel; %mm from bregma. 
yaxis_stereo = ((1:size(image_stereo,1)) - bregmaypix_stereo) * MmPerPixel; %mm from bregma. 

%% stereotaxic coordinates of the brain image
stereox = xaxis_stereo;
stereoy = yaxis_stereo;
[STEREOX, STEREOY] = meshgrid(stereox, stereoy);


%% interpolate in 2D
[STEREOX_i, STEREOY_i] = meshgrid(stereox_i, stereoy_i);

image_stereo_i = interp2(STEREOX, STEREOY, image_stereo, STEREOX_i, STEREOY_i, ...
    interpMethod, 0);%changed from nearest to spline 25/1/17

%% sanity check
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
% title('stereotaxic coords')


