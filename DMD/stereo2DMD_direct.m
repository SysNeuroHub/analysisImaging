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

%% x-y grid
xgridfrombregma = -4:.4:4; %[mm]
ygridfrombregma = -4:.4:3; %A>0, P<0

xgrid = 1/MmPerPixel_t * xgridfrombregma + bregma(2);
ygrid = -1/MmPerPixel_t * ygridfrombregma + bregma(1);


%% CCF contour projected to 2D (showAllenCCFBregmaLambda.m)

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

hline(ygrid, axpatch,'-','w');
vline(xgrid, axpatch,'-','w');

axpatch.Position = [0 0 1 1];
saveName = ['CCFBL_' num2str(width) 'x' num2str(height)];

screen2png(saveName, fpatch);
close(fpatch);
ctxOutlines = rgb2gray(imread([saveName '.png']));    
ctxOutlines = double(ctxOutlines/max(ctxOutlines(:)));
subplot(2,3,1); imagesc(ctxOutlines); hold on; plot(bregma(2), bregma(1),'ro'); axis equal tight; title('image in stereotaxic coords');

%% CCF contour in 2D 
MRIdir = '/home/daisuke/Documents/git/analysisImaging/MROIDMD';
oriimg = niftiread(fullfile(MRIdir, 'pattern_generation/Allen_annotation_modified.nii'));
oriimg_info = niftiinfo(fullfile(MRIdir, 'pattern_generation/Allen_annotation_modified.nii'));
mrsize_xy = [size(oriimg,3), size(oriimg,1)]; %[x y]

%% register CCF2D contour from UCL to individual brain

MmPerPixel_mr = 0.1; %[mm]
bregmak_c = allenCCFbregma()*MmPerPixel_mr;
bregmak = [bregmak_c(1) bregmak_c(3)]+0.5; %[y x]
movingPoints = [bregma(2) bregma(1); bregma(2)+1/MmPerPixel_t bregma(1)]; %[x y]
fixedPoints = [bregmak(2) bregmak(1); bregmak(2)+1/MmPerPixel_mr bregmak(1)]; %[x y]
tform3 = fitgeotform2d(movingPoints, fixedPoints, 'similarity');
%  [movingPoints,fixedPoints] = cpselect(ctxOutlines,ctxOutlines_DMD,'Wait',true);
%  tform3 = fitgeotrans(movingPoints,fixedPoints, 'similarity');

refDMD = imref2d(mrsize_xy);
usFactor = 5;
refDMDbig = imref2d(usFactor * mrsize_xy);
refDMDbig.XWorldLimits = refDMD.XWorldLimits;
refDMDbig.YWorldLimits = refDMD.YWorldLimits;

ctxOutlines_reg = imwarp(ctxOutlines,tform3,'linear','OutputView', refDMDbig, 'FillValues',0);
% subplot(2,3,2); imshowpair(ctxOutlines_DMD, ctxOutlines_reg); title('stereo image (m) registered to CCF (g)');
subplot(2,3,2); imagesc(ctxOutlines_reg); axis equal tight; title('stereo image registered to upscaled CCF');

%% texture mapping from 2D to 3D in Allen CCF (Kim) space
surviveR=[18 24 44 51 65 72 79 86 93 100 122 136 143 150 164 171 178 185 192 199 206 213 226 238 298 325 332 346 353 360];
surviveL=surviveR+2000; %from project_anno.m
V = flip(ismember(oriimg, [surviveL surviveR]), 3);
V = imresize3(V, usFactor, 'linear');

 [surfX, surfY, surfZ] = vol2Surf(V, 50*usFactor);
 %plot3(surfX, surfY, surfZ, '.');

% 2. Initialize texture volume
TexVol = zeros(size(V));

% 3. Define projection mapping
% Normalize x,y surface coordinates to image pixel coordinates
I = ctxOutlines_reg;
if (size(V,1) ~= size(I,2)) | (size(V,3) ~= size(I,1))
    error('volume size does not match image size');
end

u = round( rescale(surfY, 1, size(I,2), "InputMin",1,"InputMax", size(V,1)) ); %CORRECT??
v = round( rescale(surfZ, 1, size(I,1), "InputMin",1,"InputMax", size(V,3)) ); %CORRECT??

% 4. Assign projected intensity to surface voxels
for k = 1:length(surfX)
    TexVol(surfY(k), surfX(k), surfZ(k)) = I(v(k), u(k));
end

TexVol = flip(TexVol, 3); %why is this needed?
subplot(2,3,3); isosurface(TexVol); axis equal tight;
ax = gca; set(gca, 'view', [120 20]);
title('Projection mapped to upscaled CCF');

% 5. Smooth or fill small holes
TexVolSmooth = imdilate(TexVol, strel('sphere', 1));
info = oriimg_info;
info.Datatype = 'double';
info.Transform.T(1:3,1:3) = info.Transform.T(1:3,1:3)/usFactor;
info.PixelDimensions = info.PixelDimensions/usFactor;
info.ImageSize = size(TexVolSmooth);
niftiwrite(TexVolSmooth,'TexVolSmooth.nii', info);


%% warp to individual brain (in a shell script)
subjectName = 'tmpD';

copyfile(fullfile(MRIdir, 'pattern_generation/Brain_template.nii'), ...
    fullfile(MRIdir, subjectName,'Brain_template.nii'));


niftiwrite_us(fullfile(MRIdir,subjectName,'T2w_brain.nii'), usFactor);
cmdStr = [fullfile(MRIdir,'pattern_generation/AtlasTexVol_to_T2.sh') ' ' 'TexVolSmooth.nii' ' ' fullfile(MRIdir, subjectName)];
system(cmdStr); %output TexVol_T2.nii


%% project back from 3D to 2D
TexVol_T2 = niftiread('TexVolSmooth_T2.nii');

%% volume mask returns less noisier image than surface mask
%brainMask = niftiread(fullfile(MRIdir,subjectName,'T2w_brain_us.nii'));
TexImg = getSufraceData2(TexVol_T2, brainMask>0);
%[~, TexImg] = getSurfaceData3(TexVol_T2, 'last', 0);  %slightly faster but noisier
subplot(2,3,4); imagesc(TexImg);axis equal tight; title('stereo warped to T2');

%% convert to OI
load('/home/daisuke/Documents/git/analysisImaging/MROIDMD/tmpD/Atlas_reg_info.mat',...
    'tform','tform2','mrwarpedtoDMD');

% from CCF to OI
OIsize = [1080 1080];  %[y x]
DMDSize = [500 800]; %[y x]
fixedRef  = imref2d(OIsize, 0.0104, 0.0104);  % example pixel sizes in mm
movingRef = imref2d(size(TexImg), MmPerPixel_mr/usFactor,  MmPerPixel_mr/usFactor);
TexImgwarped = imwarp(TexImg,movingRef,tform,'linear','OutputView',fixedRef);
TexImgwarpedtoDMD = imwarp(TexImgwarped,tform2,'linear','OutputView',imref2d(DMDSize));

subplot(2,3,5); imagesc(TexImgwarped);axis equal tight; title('T2 warped to OI');
subplot(2,3,6); imshowpair(mrwarpedtoDMD, TexImgwarpedtoDMD);axis equal tight; title('OI warped to DMD(m)');

