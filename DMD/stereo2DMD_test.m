%% prepare DMD image

scale = 5/9;%
width = scale*1440;
height = scale*900;
brainImage = zeros(height,width);
MmPerPixel_t = 0.0104 / scale; %measured w scale 27/1/25 from getMmPerPix.m
 
xfrombregma = 3.5;%-2.7; %[mm]
yfrombregma = -3.6; %A>0, P<0
radiusmm = MmPerPixel_t*[1.5 5 10]; %[mm];
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
ctxOutlines = ctxOutlines/max(ctxOutlines(:));


%% register CCF2D contour from UCL to individual brain
load('/home/daisuke/Documents/MRI/output/tmpC/Atlas_reg_info.mat', 'mrwarpedtoDMD','TotalBrainImage');
ctxOutlines_DMD = (TotalBrainImage==0)-(mrwarpedtoDMD==0);

 tform3 = imregcorr(ctxOutlines, ctxOutlines_DMD); %NG

 % [optimizer,metric] = imregconfig("monomodal");
% tform3= imregtform(ctxOutlines, ctxOutlines_DMD,"affine",optimizer,
% metric); %NG

ctxOutlines_reg = imwarp(ctxOutlines,tform3,'cubic','OutputView',imref2d(size(ctxOutlines_DMD)));
imshowpair(ctxOutlines_DMD, ctxOutlines_reg);
