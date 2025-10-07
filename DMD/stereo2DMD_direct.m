%% prepare DMD image

scale = 5/9;%
width = scale*1440;
height = scale*900;
brainImage = zeros(height,width);
MmPerPixel_t = 0.0104 / scale; %measured w scale 27/1/25 from getMmPerPix.m
 
xfrombregma = 3.5;%-2.7; %[mm]
yfrombregma = -3.6; %A>0, P<0
bregma = [scale*(380-20) width/2+1]; 
lambda = [scale*(825-20) width/2+1];

xfrombregmapix = 1/MmPerPixel_t * xfrombregma + bregma(2);
yfrombregmapix = -1/MmPerPixel_t * yfrombregma + bregma(1);

%% CCF contour projected to 2D (showAllenCCFBregmaLambda.m)

fpatch=figure;
fpatch.InnerPosition = [1 1 width height];

image(zeros(size(brainImage)));colormap(gray);
addAllenCtxOutlines(bregma, lambda, 'w', MmPerPixel_t);%this looks at lambda and shrinks the CCF

hold on;
axis ij image off
xlim([1 size(brainImage,2)]);
ylim([1 size(brainImage,1)]);

axpatch = gca;
axpatch.Position = [0 0 1 1];
saveName = ['CCFBL_' num2str(width) 'x' num2str(height)];

screen2png(saveName, fpatch);
close(fpatch);
ctxOutlines = rgb2gray(imread([saveName '.png']));    
ctxOutlines = double(ctxOutlines/max(ctxOutlines(:)));
subplot(2,3,1); imagesc(ctxOutlines); axis equal tight; title('image in stereotaxic coords');

%% CCF contour in 2D 
MRIdir = '/home/daisuke/Documents/git/analysisImaging/MROIDMD';
oriimg = niftiread(fullfile(MRIdir, 'pattern_generation/Allen_annotation_modified.nii'));
oriimg_info = niftiinfo(fullfile(MRIdir, 'pattern_generation/Allen_annotation_modified.nii'));
[proj_anno, proj_anno_cortex, ROI_info] = project_anno(oriimg);

proj_brain=proj_anno_cortex;
Brainimage={};
for i=1:size(ROI_info,1)
    ROI_all=zeros(size(proj_brain));
    ROI_all(find(proj_brain==ROI_info{i,1}))=i;
    ROI_all2=imerode(ROI_all,strel('disk',1));
    %ROI_all3=imerode(ROI_all2,strel('diamond',1));
    Brainimage{i}=ROI_all2;
end

ImageRow=[];
ROI_information=ROI_info;
for re=1:size(Brainimage,2)
    ImageRow=cat(3,Brainimage{re},ImageRow);
    ROI_information{re,1}=re;
end
TotalBrainImage=sum(ImageRow,3);

ctxOutlines_DMD = (TotalBrainImage==0)-(proj_brain==0);


%% register CCF2D contour from UCL to individual brain

 % tform3 = imregcorr(ctxOutlines, ctxOutlines_DMD); %NG

%  [optimizer,metric] = imregconfig("multimodal");
%  fixedRef  = imref2d(size(ctxOutlines_DMD), 0.1, 0.1);  % example pixel sizes in mm
%  movingRef = imref2d(size(ctxOutlines), MmPerPixel_t, MmPerPixel_t);
% 
% tform3= imregtform(ctxOutlines, movingRef, ctxOutlines_DMD, fixedRef,"similarity",optimizer, metric); %NG
% tform3= imregtform(imfill(ctxOutlines, 'holes'), movingRef, imfill(ctxOutlines_DMD,'holes'), fixedRef,"similarity",optimizer, metric); %NG

  [movingPoints,fixedPoints] = cpselect(ctxOutlines,ctxOutlines_DMD,'Wait',true);
  tform3 = fitgeotrans(movingPoints,fixedPoints, 'similarity');

refDMD = imref2d(size(ctxOutlines_DMD));
usFactor = 5;
refDMDbig = imref2d(usFactor * size(ctxOutlines_DMD));
refDMDbig.XWorldLimits = refDMD.XWorldLimits;
refDMDbig.YWorldLimits = refDMD.YWorldLimits;

%ctxOutlines_reg = imwarp(ctxOutlines,tform3,'cubic','OutputView', refDMD);
ctxOutlines_reg = imwarp(ctxOutlines,tform3,'linear','OutputView', refDMDbig, 'FillValues',0);
% subplot(2,3,2); imshowpair(ctxOutlines_DMD, ctxOutlines_reg); title('stereo image (m) registered to CCF (g)');
subplot(2,3,2); imagesc(ctxOutlines_reg); axis equal tight; title('stereo image registered to upscaled CCF');

%% texture mapping from 2D to 3D in Allen CCF (Kim) space
%V = flip((oriimg~=0).*(oriimg~=2000), 3);
surviveR=[18 24 44 51 65 72 79 86 93 100 122 136 143 150 164 171 178 185 192 199 206 213 226 238 298 325 332 346 353 360];
surviveL=surviveR+2000; %from project_anno.m
V = flip(ismember(oriimg, [surviveL surviveR]), 3);
%V = imfill(V);
V = imresize3(V, usFactor, 'linear');

 [surfX, surfY, surfZ] = vol2Surf(V, 50*usFactor);

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
title('Stereo image projection mapped to CCF');

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

% if ~isempty(angle)
%     TexVol_T2 = imrotate3(TexVol_T2, angle(1), [1 0 0],'linear','crop'); %roll
%     TexVol_T2 = imrotate3(TexVol_T2, angle(2), [0 1 0],'linear','crop'); %pitch
%     TexVol_T2 = imrotate3(TexVol_T2, angle(3), [0 0 1],'linear','crop'); %yaw
% end
%
%[~, TexImg] = getSurfaceData(TexVol_T2); % extremely noisy

%% volume mask returns less noisier image than surface mask
brainMask = niftiread(fullfile(MRIdir,subjectName,'T2w_brain_us.nii'));
TexImg = getSufraceData2(TexVol_T2, brainMask>0);
subplot(2,3,4); imagesc(TexImg);axis equal tight; title('stereo warped to T2');

%% convert to OI
load('/home/daisuke/Documents/git/analysisImaging/MROIDMD/tmpD/Atlas_reg_info.mat',...
    'tform','tform2','mrwarpedtoDMD');

% from CCF to OI
OIsize = [1080 1080];
fixedRef  = imref2d(OIsize, 0.0104, 0.0104);  % example pixel sizes in mm
movingRef = imref2d(size(TexImg), 0.1/usFactor,  0.1/usFactor);
TexImgwarped = imwarp(TexImg,movingRef,tform,'linear','OutputView',fixedRef);
TexImgwarpedtoDMD = imwarp(TexImgwarped,tform2,'linear','OutputView',imref2d([500 800]));

subplot(2,3,5); imagesc(TexImgwarped);axis equal tight; title('T2 warped to OI');
subplot(2,3,6); imshowpair(mrwarpedtoDMD, TexImgwarpedtoDMD);axis equal tight; title('OI warped to DMD(m)');

