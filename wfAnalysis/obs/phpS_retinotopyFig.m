addpath(genpath('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\dsbox'));

prefMap = imresize(prefMap,size(mimg));

%prepare brain image ranging 0-1
maxprctile = 90;
minprctile = 25;
tgtimg = log(mimg);
brainimage = (tgtimg - prctile(tgtimg(:),minprctile))...
    /(prctile(tgtimg(:),maxprctile) - prctile(tgtimg(:),minprctile));
brainimage(brainimage>1)=1;
brainimage(brainimage<0)=0;
imagesc(brainimage);colorbar;


imagesc(prefMap,'alphadata',brainimage);
set(gca,'color','k');
caxis([1 p.nstim-1]);
cb=colorbar;
cb.Ticks = 1:p.nstim-1;
cb.TickLabels = stimSequence.labels(cb.Ticks);
colormap(jet);
axis equal tight;

%scale bar
pixsize = 5e-3;%pixel size of the camera sensor in mm
binning = 2; %binning factor 
magnification = 0.5;%objective
npixpermm = magnification/binning/pixsize;
line([5 5+npixpermm],[5 5],'color','w','linewidth',2);

savePaperFigure(gcf,'phpS-2020-07-06_2-1','k');
