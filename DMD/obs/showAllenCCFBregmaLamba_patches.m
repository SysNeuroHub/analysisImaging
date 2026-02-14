addpath(genpath('C:\Users\dshi0006\allenCCF'));
addpath(genpath('C:\Users\dshi0006\git\analysisImaging'));
addpath('C:\Documents\git\dsbox\visualization');

scale = 1/3;

width =  400;
height = 300;
bregma = [scale*(380-20) width/2+1]+0.5; %[y x]
lambda = [scale*(825-20) width/2+1]+0.5;

brainImage = zeros(height,width);
% MmPerPixel_t = 6.5e-3/0.5; %nominal value
%MmPerPixel_t = 0.75*6.5e-3/0.5; %educated guess
MmPerPixel_img = 0.0104 / scale; %measured w scale 27/1/25 from getMmPerPix.m

f=figure;
f.InnerPosition = [1 1 width height];
%rectangle('position',f.InnerPosition,'edgecolor','r');hold on; %otherwise black region outside brain will be trimmed in savePaperFigure
image(zeros(size(brainImage)));colormap(gray);

addAllenCtxOutlines(bregma, lambda, 'w', MmPerPixel_img);%this looks at lambda and shrinks the CCF
scatter(bregma(2), bregma(1),600,'markerEdgecolor','w','MarkerFaceColor','w','marker','x');
scatter(lambda(2), lambda(1),600,'markerEdgecolor','w','MarkerFaceColor','w','marker','x');
axis ij image off
xlim([1 size(brainImage,2)]);
ylim([1 size(brainImage,1)]);
line([width/2-50 width/2-50+1/MmPerPixel_img], [800 800],'linewidth',2,'color','w');

ax = gca;
ax.Position = [0 0 1 1];

%% superimpose grid on CCF
xfrombregma = -4:1:4; %[mm]
yfrombregma = -4:1:3; %A>0, P<0

xgrid = 1/MmPerPixel_img * xfrombregma + bregma(2);
ygrid = -1/MmPerPixel_img * yfrombregma + bregma(1);

hline(ygrid, ax,'-','w');
vline(xgrid, ax,'-','w');


saveName = ['CCFBL_' num2str(width) 'x' num2str(height) 'pix_' num2str(numel(xgrid)-1) 'x' num2str(numel(ygrid)-1) 'grid'];
saveServer = '~/tmp';
saveDir = fullfile(saveServer, saveName);
mkdir(saveDir);
mkdir(fullfile(saveDir, 'stereo'));
screen2png(fullfile(saveDir, [saveName '_stereo']), f);
close(f);

%% show only patches along the grid

gridNumber = 0;
imageStereo = [];
imageName = cell(1);
for yy = 1:numel(yfrombregma)-1
    for xx = 1:numel(xfrombregma)-1
        gridNumber = gridNumber + 1;

        fpatch=figure;
        fpatch.InnerPosition = [1 1 width height];

        image(zeros(size(brainImage)));colormap(gray);
        addAllenCtxOutlines(bregma, lambda, 'w', MmPerPixel_img);%this looks at lambda and shrinks the CCF
        hold on;
        rectangle('position',[xgrid(xx) ygrid(yy+1) ...
            diff(xgrid(xx:xx+1)) -diff(ygrid(yy:yy+1))],...
            'facecolor','w');
        axis ij image off
        xlim([1 size(brainImage,2)]);
        ylim([1 size(brainImage,1)]);
        
        axpatch = gca;
        axpatch.Position = [0 0 1 1];

        thisName = [num2str(gridNumber) '_wCCF'];
        screen2png(fullfile(saveDir, 'stereo', thisName),fpatch);
        close(fpatch);

        
        %% load
        imageLoaded = rgb2gray(imread(fullfile(saveDir, 'stereo', [thisName '.png'])));    
        imageLoaded = double(imageLoaded/max(imageLoaded(:)));

        imageStereo = cat(3, imageStereo, imageLoaded);
        imageName{gridNumber} = thisName; 
    end
end

%% save everything in one file
save(fullfile(saveDir, [saveName '_stereo']), 'imageStereo','bregma','lambda', 'MmPerPixel_img');



%% show only CCF without patches
fpatch=figure;
fpatch.InnerPosition = [1 1 width height];

image(zeros(size(brainImage)));colormap(gray);
addAllenCtxOutlines(bregma, lambda, 'w', MmPerPixel_img);%this looks at lambda and shrinks the CCF

axis ij image off
xlim([1 size(brainImage,2)]);
ylim([1 size(brainImage,1)]);

axpatch = gca;
axpatch.Position = [0 0 1 1];

thisName = ['CCF'];
screen2png(fullfile(saveDir, 'stereo', thisName),fpatch);
close(fpatch);


imageLoaded = rgb2gray(imread(fullfile(saveDir, 'stereo', [thisName '.png'])));
imageLoaded = double(imageLoaded/max(imageLoaded(:)));

imageStereo_CCF = imageLoaded;

save(fullfile(saveDir, [saveName '_stereo']), 'imageStereo_CCF','-append');
