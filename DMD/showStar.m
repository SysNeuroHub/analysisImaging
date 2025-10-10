addpath(genpath('C:\Users\dshi0006\git\analysisImaging'));

scale = 5/9;%FIXED
width = scale*1440;
height = scale*900;
brainImage = zeros(height,width);
MmPerPixel_t = 0.0104 / scale; %measured w scale 27/1/25 from getMmPerPix.m

f=figure;
f.InnerPosition = [1 1 width height];
image(zeros(size(brainImage)));colormap(gray); 
hold on;
line([width/2-scale*50 width/2-scale*50+1/MmPerPixel_t], scale*[800 800],'linewidth',2,'color','w');

%% add reference star
ax = gca;
mark=plot(width/2+1, height/2+1,'wp','MarkerSize',0.6*height,'MarkerFaceColor','none','MarkerEdgeColor','w');
vline(width/2+1, ax, '-','w');
hline(height/2+1, ax, '-','w');
text(width/4, height/4,'TL','FontSize',30,'color','w')
ax.Position = [0 0 1 1];

saveName = ['star_' num2str(width) 'x' num2str(height)];
screen2png(saveName,f);%
