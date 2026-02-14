% run this function and take two images when
% OASIS macro setup has been modified (polygon orientation, optical filters and dichroics)
% created from saveRefImg.m

% addpath(genpath('C:\Users\dshi0006\allenCCF'));
% addpath(genpath('C:\Users\dshi0006\git\analysisImaging'));

gitDir = '/home/daisuke/Documents/git/analysisImaging';
refDir = fullfile(gitDir, '/DMD/references');
addpath(genpath(gitDir));

refdate = '20260214';

%% take two images with camera, full ROI, without binning. bottom dichro 465, top dichro off
%1) projection zone: all DMD mirrors ON. 
% save as fullfile(gitDir, refDir, 'DMDprojectionZone', ['DMDprojectionZone' date  '.tif'])
%
%2) reference star, project star_800x500.png, created by fullfile(gitDir,'DMD',' showStar')
% save as fullfile(gitDir, 'MROIDMD/matlab_functions', ['star_' num2str(camImg.imageSize(2)) 'x' num2str(camImg.imageSize(1)) '.tif']) (OBSOLETE)
% save as fullfile(gitDir, refDir, 'star', ['star'  date '.tif'])


%% camera image info
binning = 2; % scale = 1/binning
width = 1168/binning;
height = 900/binning;
bregma = [360/binning width/binning+1];
lambda = [805/binning width/binning+1];
MmPerPixel = 0.0098 * binning; %measured w scale 14/2/26 from getMmPerPix.m
topLeftPix = [288 1+8*25]; %column value minimum difference is 8
camImg = cameraImageInfo([height width], bregma, lambda, MmPerPixel,binning, topLeftPix);

save(fullfile(refDir, ['camImg_' refdate]),'camImg');


%% save CCF image with bregma, lambda, scale and projection zone
saveName = ['CCFBL_' num2str(camImg.imageSize(2)) 'x' num2str(camImg.imageSize(1)) 'pix_top' num2str(topLeftPix(1)) '_left' num2str(topLeftPix(2))];
saveDir = fullfile(refDir, 'DMD images', saveName);
mkdir(saveDir);

f = showRefImg(camImg, fullfile(refDir, 'DMDprojectionZone'), ['DMDprojectionZone' refdate  '.tif']);
exportPng4DMD(fullfile(saveDir, [saveName '_' refdate]), f, 1);
close(f);
