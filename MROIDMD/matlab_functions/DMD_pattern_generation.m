function ImageOut = DMD_pattern_generation(ROI,B,C,proj_brain,ROI_info, mapconfwarpedtoDMD)
% Input: (ROI,B,C, proj_brain, ROI_info, mapconfwarpedtoDMD)
%       ROI: [row vector of ROIs] [1 60]
%       B: Gaussian smooting factor (sd), default=10
%       C: erosion factor (disk & diamond), default=1
%       proj_brain, ROI_info, mapconfwarpedtoDMD: output of MR2DMD, saved in 'Atlas_reg_info.mat'
%
% created from  matlab_functions/DMD_pattern_generation

showFig = false;

if(~exist('B','var'))
    B=10;
    fprintf('B not inserted; Gaussian smoothing SD=10 \n')
end

if(~exist('C','var'))
    C=1;
    fprintf('C not inserted; erosion factor=1\n')
end

if nargin < 6
    mapconfwarpedtoDMD = [];
end


%% Brainimage: DMD image 500x800pix of all CCF regions
    Brainimage={};
    for i=1:size(ROI_info,1)
        ROI_all=zeros(size(proj_brain));
        ROI_all(find(proj_brain==ROI_info{i,1}))=1;
        Brainimage{i}=ROI_all;
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

if showFig
    jetkey=jet(400);
    figure;subplot(2,1,1);imagesc(ImageTotal);hold on; 
   if ~isempty(mapconfwarpedtoDMD)
       contour(mapconfwarpedtoDMD,'r');
   end
    title('all CCF regions');

    subplot(2,1,2);imagesc(2*Imageraw+ImageOut);colormap([1 1 1;jetkey])
    titlename=sprintf('ROI%d: (%d/2560*1600 pix: %d (percent))',ROI(1),length(find(ImageOut)),100*length(find(ImageOut))/(2560*1600));
    title(titlename)
end

% Outputname=sprintf('ROI_%d.png',ROI(1));
% imwrite(ImageOut,Outputname)
% fprintf('Done! ROI image file has been created.\n')














