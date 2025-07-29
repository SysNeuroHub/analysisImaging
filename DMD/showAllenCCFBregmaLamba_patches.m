addpath(genpath('C:\Users\dshi0006\allenCCF'));
addpath(genpath('C:\Users\dshi0006\git\analysisImaging'));

%load('M:\Subjects\himiko\2025-01-23_1\dataSummary_amber.mat', 'dataSummary');
%brainImage = dataSummary.meanImage;
width = 1168;%1000;
height = 900;
brainImage = zeros(height,width);
% MmPerPixel_t = 6.5e-3/0.5; %nominal value
%MmPerPixel_t = 0.75*6.5e-3/0.5; %educated guess
MmPerPixel_t = 0.0104; %measured w scale 27/1/25 from getMmPerPix.m

f=figure;
f.InnerPosition = [1 1 width height];
%rectangle('position',f.InnerPosition,'edgecolor','r');hold on; %otherwise black region outside brain will be trimmed in savePaperFigure
image(zeros(size(brainImage)));colormap(gray);

bregma = [380 width/2+1];
lambda = [825 width/2+1];
% bregma = [380 515];%[yy xx] himiko
% lambda = [825 515];%[yy xx] himiko
addAllenCtxOutlines(bregma, lambda, 'w', MmPerPixel_t);%this looks at lambda and shrinks the CCF
scatter(bregma(2), bregma(1),600,'markerEdgecolor','w','MarkerFaceColor','w','marker','x');
scatter(lambda(2), lambda(1),600,'markerEdgecolor','w','MarkerFaceColor','w','marker','x');
axis ij image off
xlim([1 size(brainImage,2)]);
ylim([1 size(brainImage,1)]);
line([width/2-50 width/2-50+1/MmPerPixel_t], [800 800],'linewidth',2,'color','w');

ax = gca;
ax.Position = [0 0 1 1];

saveName = ['CCFBL_' num2str(width) 'x' num2str(height)];
screen2png(saveName,f);%
% savePaperFigure(f, saveName);%use .eps to paint a ROI then save as a png
%exportgraphics(f,'CCFBL_1000x900.png','BackgroundColor','k'); %does not preserve original pixel size

%% superimpose grid on CCF
xfrombregma = -4:1:0; %[mm]
yfrombregma = -4:3; %A>0, P<0

xgrid = 1/MmPerPixel_t * xfrombregma + bregma(2);
ygrid = -1/MmPerPixel_t * yfrombregma + bregma(1);

hline(ygrid, ax,'-','w');
vline(xgrid, ax,'-','w');

screen2png([saveName '_grid'],f);%


%% show only patches along the grid

gridNumber = 0;
for yy = 1:numel(yfrombregma)-1
    for xx = 1:numel(xfrombregma)-1
        gridNumber = gridNumber + 1;

        fpatch=figure;
        fpatch.InnerPosition = [1 1 width height];

        image(zeros(size(brainImage)));colormap(gray);
        addAllenCtxOutlines(bregma, lambda, 'w', MmPerPixel_t);%this looks at lambda and shrinks the CCF
        hold on;
        rectangle('position',[xgrid(xx) ygrid(yy+1) ...
            diff(xgrid(xx:xx+1)) -diff(ygrid(yy:yy+1))],...
            'facecolor','w');
        axis ij image off
        xlim([1 size(brainImage,2)]);
        ylim([1 size(brainImage,1)]);
        
        axpatch = gca;
        axpatch.Position = [0 0 1 1];
        screen2png([num2str(gridNumber) '_wCCF'],fpatch);
        close(fpatch);
    end
end