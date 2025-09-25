% Requirements: Atlas_anno_to_T2.nii
% Input: (ROI,B,C)
%       ROI: [row vector of ROIs]
%       B: Gaussian smooting factor (sd), default=10
%       C: erosion factor (disk & diamond), default=1
%       D: 1: Pop up Total Brain Image figure; else: X

function DMD_pattern_generation(ROI,B,C,D)

%% Load data
if(~exist('B','var'))
    B=10;
    fprintf('B not inserted; Gaussian smoothing SD=10 \n')
end

if(~exist('C','var'))
    C=1;
    fprintf('C not inserted; erosion factor=1\n')
end

if(~exist('D','var'))
    D=0;
    fprintf('D not inserted; if want to pop up ImageTotal, insert D\n')
end


load('Atlas_reg_info.mat')

%% Pattern generation

Brainfname=sprintf('Brainimage_%d_%d.mat',B,C);

if(~exist(Brainfname,'file'))
fprintf('Brainimage does not exist... creating new Brainimage...\n')

jetkey=jet(400);
if D==1
figure;imagesc(TotalBrainImage);colormap([1 1 1;jetkey])
end


Brainimage={};
for i=1:size(ROI_info,1)
    ROI_all=zeros(size(proj_brain));
    ROI_all(find(proj_brain==ROI_info{i,1}))=1;
    Brainimage{i}=ROI_all;
end

fprintf('Brainimage is created... Saving Brainimage...\n')
save(Brainfname,'Brainimage','-v7.3')

else
fprintf('Brainimage exists...\n')
load(Brainfname)
end

ImageRow=zeros(size(Brainimage{1}));
for re=1:size(Brainimage,2)
    ImageRow2=imerode(Brainimage{re},strel('diamond',2));
    ImageRow=ImageRow+re*ImageRow2;
end
ImageTotal=ImageRow;

ImageRow=zeros(size(Brainimage{1}));
for re=ROI
    ImageRow=ImageRow+Brainimage{re};
end
gaussimg=imgaussfilt(ImageRow,B);
gaussimg(gaussimg<0.7)=0;
gaussimg(gaussimg>0)=1;

ImageRow2=imerode(gaussimg,strel('disk',C));
ImageRow3=imerode(gaussimg,strel('diamond',C)); 
ImageOut=ImageRow3;

Imageraw=zeros(size(Brainimage{1}));
for i = 1:size(Brainimage,2)
    Imagebr = i*Brainimage{i};
    Imageraw = Imageraw+Imagebr;
end
Imageraw=Imageraw-imerode(Imageraw,strel('diamond',2));
Imageraw(Imageraw>0)=1;

jetkey=jet(400);
figure;subplot(2,1,1);imagesc(ImageTotal+mapconfwarpedtoDMD)
subplot(2,1,2);imagesc(2*Imageraw+ImageOut);colormap([1 1 1;jetkey])
titlename=sprintf('ROI%d: (%d/2560*1600 pix: %d (percent))',ROI(1),length(find(ImageOut)),100*length(find(ImageOut))/(2560*1600));
title(titlename)

Outputname=sprintf('ROI_%d.png',ROI(1));
imwrite(ImageOut,Outputname)
fprintf('Done! ROI image file has been created.\n')














