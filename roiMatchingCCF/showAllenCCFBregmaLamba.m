setPath_analysisImaging;

expt.subject = 'saruta';
expt.expDate = '2025-03-14_2';
expt.expNum = 1;

load(fullfile(rawDataDir, expt.subject, expt.expDate, num2str(expt.expNum), 'dataSummary_amber.mat'));
brainImage = dataSummary.meanImage;
width = size(brainImage,2);%1000;
height = size(brainImage,1);%900;
MmPerPixel_t = 0.0104; %measured w scale 27/1/25 from getMmPerPix.m

f=figure;
f.InnerPosition = [1 1 width height];
imagesc(brainImage);colormap(gray);

hold on;
bregma = [];%[377 507];
lambda = [];%[838 510];
addAllenCtxOutlines(bregma, lambda, 'w', MmPerPixel_t);%this looks at lambda and shrinks the CCF
axis image off
xlim([1 size(brainImage,2)]);
ylim([1 size(brainImage,1)]);
margin = 50;%pix
line([margin margin+1/MmPerPixel_t], [margin margin],'linewidth',2,'color','w');

ax = gca;
ax.Position = [0 0 1 1];

saveName = [expt.subject '_CCFBL_' num2str(width) 'x' num2str(height)];
screen2png(saveName,f);

