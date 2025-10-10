%% prepare DMD image

scale = 5/9;%
width = scale*1440;
height = scale*900;
brainImage = zeros(height,width);
MmPerPixel_t = 0.0104 / scale; %measured w scale 27/1/25 from getMmPerPix.m
 
xfrombregma = 3.5;%-2.7; %[mm]
yfrombregma = -3.6; %A>0, P<0
bregma = [scale*(380-20) width/2+1]+0.5; %[y x]
lambda = [scale*(825-20) width/2+1]+0.5;

% xfrombregmapix = 1/MmPerPixel_t * xfrombregma + bregma(2);
% yfrombregmapix = -1/MmPerPixel_t * yfrombregma + bregma(1);

% x-y grid
xgridfrombregma = -4:.4:4; %[mm]
ygridfrombregma = -4:.4:3; %A>0, P<0

xgrid = 1/MmPerPixel_t * xgridfrombregma + bregma(2);
ygrid = -1/MmPerPixel_t * ygridfrombregma + bregma(1);



fpatch=figure;
fpatch.InnerPosition = [1 1 width height];

image(zeros(size(brainImage)));colormap(gray);
addAllenCtxOutlines(bregma, lambda, 'w', MmPerPixel_t);%this looks at lambda and shrinks the CCF

hold on;

scatter(bregma(2), bregma(1), scale*600,'markerEdgecolor','w','MarkerFaceColor','w','marker','x');

axis ij image off
xlim([1 size(brainImage,2)]);
ylim([1 size(brainImage,1)]);

axpatch = gca;
% hline(ygrid, axpatch,'-','w');
% vline(xgrid, axpatch,'-','w');

axpatch.Position = [0 0 1 1];
saveName = ['CCFBL_' num2str(width) 'x' num2str(height)];

screen2png(saveName, fpatch);
close(fpatch);
ctxOutlines = rgb2gray(imread([saveName '.png']));    
ctxOutlines = double(ctxOutlines/max(ctxOutlines(:)));
subplot(2,3,1); imagesc(ctxOutlines); hold on; plot(bregma(2), bregma(1),'ro'); axis equal tight; title('image in stereotaxic coords');

%% CCF contour in 2D 
if ispc
MRIdir = 'C:/Documents/git/analysisImaging/MROIDMD';
else
 MRIdir = '/home/daisuke/Documents/git/analysisImaging/MROIDMD';
end
% V = flip(niftiread(fullfile(MRIdir, 'pattern_generation/Brain_template.nii'))>0,3);
V = niftiread(fullfile(MRIdir, 'pattern_generation/Brain_template.nii'))>0;
V_info = niftiinfo(fullfile(MRIdir, 'pattern_generation/Brain_template.nii'));
mrsize_xy = [size(V,3), size(V,1)]; %[x y]

%% register CCF2D contour from UCL to individual brain

MmPerPixel_mr = 0.1; %[mm]
bregmak_c = allenCCFbregma()*MmPerPixel_mr;
bregmak = [bregmak_c(1) bregmak_c(3)]+0.5; %[y x]
movingPoints = [bregma(2) bregma(1); bregma(2)+1/MmPerPixel_t bregma(1)]; %[x y]
fixedPoints = [bregmak(2) bregmak(1); bregmak(2)+1/MmPerPixel_mr bregmak(1)]; %[x y]
tform3 = fitgeotform2d(movingPoints, fixedPoints, 'similarity'); %once

refDMD = imref2d(mrsize_xy);
usFactor = 5;
refDMDbig = imref2d(usFactor * mrsize_xy);
refDMDbig.XWorldLimits = refDMD.XWorldLimits;
refDMDbig.YWorldLimits = refDMD.YWorldLimits;

ctxOutlines_reg = imwarp(ctxOutlines,tform3,'linear', ...
    'OutputView', refDMDbig, 'FillValues',0); %all
% subplot(2,3,2); imshowpair(ctxOutlines_DMD, ctxOutlines_reg); title('stereo image (m) registered to CCF (g)');
subplot(2,3,2); imagesc(ctxOutlines_reg); axis equal tight; title('stereo image registered to upscaled CCF');

%% texture mapping from 2D to 3D in Allen CCF (Kim) space
V = imresize3(V, usFactor, 'linear');%once

[surfDepth] = getSurfaceData3(V, 'last', 0);%once
TexVol = paintSurfaceToVolume(surfDepth, ctxOutlines_reg, size(V));%all

subplot(2,3,3); isosurface(TexVol); axis equal tight;
 % f = gcf;
 % f.CurrentAxes.ZDir = 'Reverse'; %should NOT reverse
ax = gca; set(gca, 'view', [120 20]);
title('Projection mapped to upscaled CCF');

% % 5. Smooth or fill small holes?
% TexVolSmooth = imdilate(TexVol, strel('sphere', 1));
info = V_info;
info.Datatype = 'double';
info.Transform.T(1:3,1:3) = info.Transform.T(1:3,1:3)/usFactor;
info.PixelDimensions = info.PixelDimensions/usFactor;
info.ImageSize = size(TexVol);
% niftiwrite(TexVolSmooth,'TexVolSmooth.nii', info);
niftiwrite(TexVol,'TexVol.nii', info); %all


%% warp to individual brain (in a shell script)
subjectName = 'tmpD';

copyfile(fullfile(MRIdir, 'pattern_generation/Brain_template.nii'), ...
    fullfile(MRIdir, subjectName,'Brain_template.nii'));


niftiwrite_us(fullfile(MRIdir,subjectName,'T2w_brain.nii'), usFactor);

cmdStr = [fullfile('/home/daisuke/Documents/git/analysisImaging/DMD', '/AtlasTexVol_to_T2.sh') ' ' 'TexVol.nii' ' ' fullfile(MRIdir, subjectName)];
system(cmdStr); %output TexVol_T2.nii %all


%% project back from 3D to 2D
loadDir = '/home/daisuke/Documents/git/analysisImaging/MROIDMD/tmpD/';
load(fullfile(loadDir, 'Atlas_reg_info.mat'),...
    'tform','tform2','mrwarpedtoDMD', 'mrwarped','image2','mrangle');

TexVol_T2 = niftiread('TexVol_T2.nii'); %this looks correct

if ~isempty(mrangle)
    TexVol_T2 = imrotate3(TexVol_T2, mrangle(1), [1 0 0],'linear','crop'); %roll
    TexVol_T2 = imrotate3(TexVol_T2, mrangle(2), [0 1 0],'linear','crop'); %pitch
    TexVol_T2 = imrotate3(TexVol_T2, mrangle(3), [0 0 1],'linear','crop'); %yaw
end

%% volume mask returns less noisier image than surface mask
% brainMask = niftiread(fullfile(MRIdir,subjectName,'T2w_brain_us.nii'));
% TexImg = getSufraceData2(TexVol_T2, brainMask>0);
[~, TexImg] = getSurfaceData3(TexVol_T2, 'last', 0);  %slightly faster but noisier

%% convert to OI
mr_brain = niftiread(fullfile(loadDir,  'T2w_brain.nii'));
if ~isempty(mrangle)
    mr_brain = imrotate3(mr_brain, mrangle(1), [1 0 0],'linear','crop'); %roll
    mr_brain = imrotate3(mr_brain, mrangle(2), [0 1 0],'linear','crop'); %pitch
    mr_brain = imrotate3(mr_brain, mrangle(3), [0 0 1],'linear','crop'); %yaw
end
surfDepth = getSurfaceData3(mr_brain, 'last');

% subplot(2,3,4); imagesc(TexImg);axis equal tight; title('stereo warped to T2');
% volumeViewer((imresize3(mr_brain,5)>0)+10*TexVol_T2) this looks ok
subplot(2,3,4); imshowpair(imresize(surfDepth,usFactor), TexImg);axis equal tight; title('stereo warped to T2'); %OK
% but this is not correct

% from CCF to OI
OIsize = size(image2);  %[y x]
DMDSize = [500 800]; %[y x] FIXED
fixedRef  = imref2d(OIsize, 0.0104, 0.0104);  % example pixel sizes in mm
movingRef = imref2d(size(TexImg), MmPerPixel_mr/usFactor,  MmPerPixel_mr/usFactor);
TexImgwarped = imwarp(TexImg,movingRef,tform,'linear','OutputView',fixedRef);
TexImgwarpedtoDMD = imwarp(TexImgwarped,tform2,'linear','OutputView',imref2d(DMDSize));

subplot(2,3,5); imshowpair(mrwarped, TexImgwarped);axis equal tight; title('T2 warped to OI'); %NG
subplot(2,3,6); imshowpair(mrwarpedtoDMD, TexImgwarpedtoDMD);axis equal tight; title('OI warped to DMD(m)');

