addpath(genpath('C:\Users\dshi0006\allenCCF'));
% load('M:\Subjects\himiko\2025-01-23_1\dataSummary_amber.mat', 'dataSummary');
load('M:\Subjects\izanami\2025-01-30_2\dataSummary_amber.mat', 'dataSummary')
brainImage = dataSummary.meanImage;
MmPerPixel_t = 6.5e-3/0.5;

% t = Tiff('\\vault-v2.erc.monash.edu.au\MNHS-dshi0006\Subjects\L4GCaMP6s_250\2020-09-21\blue nocone.tif','r');
% blue = read(t);
% blue_t = blue;
% blue_t = fliplr(blue_t');
% blue_t = double(blue_t);
% 
% StereotaxicInfo_t.yawAngle = -5; 
% StereotaxicInfo_t.rollAngle= 20;
% StereotaxicInfo_t.pitchAngle = 0;
% StereotaxicInfo_t.bregmaxpix = 485;
% StereotaxicInfo_t.bregmaypix = 315;
% MmPerPixel_t = 1e-3*5/0.5*2;
% 
% %% option 1: show as stereotaxic coords (Shimaoka Cell Rep / elife)
% stereox_i = -8:.01:1;%-5:.01:1;
% stereoy_i = -3:.01:8;%-3:.01:6;
% blue_t_s = impix2stereo_test(blue_t, StereotaxicInfo_t, MmPerPixel_t, ...
%     stereox_i, stereoy_i);
% subplot(121);
% imagesc(stereox_i, stereoy_i, blue_t_s);
% axis equal tight
% drawTopDownCtx;
% title('show as stereotaxic coords');

%% option 2: show as pixel coord 
[blue_s_t, xaxis_stereo, yaxis_stereo] = impix2stereo_test2(brainImage, StereotaxicInfo_t, MmPerPixel_t);
subplot(122);
imagesc(blue_s_t);
axis equal tight
grid on

hold on;
bregma = [];
lambda = [];
addAllenCtxOutlines(bregma, lambda, 'w', MmPerPixel_t);%this looks at lambda and shrinks the CCF
title('show as pixels coords');

screen2png('addAllenCtxOutlines_test');

