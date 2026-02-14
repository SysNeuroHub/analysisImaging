addpath(genpath('C:\Users\dshi0006\allenCCF'));
addpath(genpath('C:\Users\dshi0006\git\analysisImaging'));

%% camera image
binning = 2; % scale = 1/binning
width = 1168/binning;
height = 900/binning;
bregma = [360/binning width/binning+1];
lambda = [805/binning width/binning+1];
MmPerPixel = 0.0104 * binning; %measured w scale 27/1/25 from getMmPerPix.m
%topLeftPix = [270 186]; %9/2/26
topLeftPix = [288 1+8*25]; %column value minimum difference is 8
camImg = cameraImageInfo([height width], bregma, lambda, MmPerPixel,binning, topLeftPix);


%% save CCF image with bregma, lambda and scale
saveName = ['CCFBL_' num2str(camImg.imageSize(2)) 'x' num2str(camImg.imageSize(1)) 'pix_top' num2str(topLeftPix(1)) '_left' num2str(topLeftPix(2))];
saveServer = '/mnt/dshi0006_market/DMD images/';%'~/tmp';
saveDir = fullfile(saveServer, saveName);
mkdir(saveDir);

f = showRefImg(camImg);
exportPng4DMD(fullfile(saveDir, saveName), f, 1);
close(f);
