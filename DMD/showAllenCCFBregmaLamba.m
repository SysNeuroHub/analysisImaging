addpath(genpath('C:\Users\dshi0006\allenCCF'));
addpath(genpath('C:\Users\dshi0006\git\analysisImaging'));

%load('M:\Subjects\himiko\2025-01-23_1\dataSummary_amber.mat', 'dataSummary');
%brainImage = dataSummary.meanImage;
width = 1000;
height = 900;
brainImage = zeros(height,width);
% MmPerPixel_t = 6.5e-3/0.5; %nominal value
%MmPerPixel_t = 0.75*6.5e-3/0.5; %educated guess
MmPerPixel_t = 0.0104; %measured w scale 27/1/25 from getMmPerPix.m

f=figure;
f.InnerPosition = [1 1 width height];
image(zeros(size(brainImage)));colormap(gray);

hold on;
bregma = [380 501];
lambda = [825 501];
% bregma = [380 515];%[yy xx] himiko
% lambda = [825 515];%[yy xx] himiko
addAllenCtxOutlines(bregma, lambda, 'w', MmPerPixel_t, 1);%this looks at lambda and shrinks the CCF
axis image off
xlim([1 size(brainImage,2)]);
ylim([1 size(brainImage,1)]);
margin = 50;%pix
line([margin margin+1/MmPerPixel_t], [margin margin],'linewidth',2,'color','w');

ax = gca;
ax.Position = [0 0 1 1];

screen2png('CCFBL_1000x900',f);
%exportgraphics(f,'CCFBL_1000x900.png','BackgroundColor','k'); %does not preserve original pixel size

