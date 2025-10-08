function [tform3, surfDepth, sizeVus] = registerStereo2CCF(bregma, MmPerPixel, usFactor)

if isempty(usFactor)
    usFactor = 5;
end

MRIdir = 'C:/Documents/git/analysisImaging/MROIDMD';
MmPerPixel_mr = 0.1; %[mm]


V = flip(niftiread(fullfile(MRIdir, 'pattern_generation/Brain_template.nii'))>0,3);
%V_info = niftiinfo(fullfile(MRIdir, 'pattern_generation/Brain_template.nii'));
mrsize_xy = [size(V,3), size(V,1)]; %[x y]


bregmak_c = allenCCFbregma()*MmPerPixel_mr;
bregmak = [bregmak_c(1) bregmak_c(3)]+0.5; %[y x]
movingPoints = [bregma(2) bregma(1); bregma(2)+1/MmPerPixel bregma(1)]; %[x y]
fixedPoints = [bregmak(2) bregmak(1); bregmak(2)+1/MmPerPixel bregmak(1)]; %[x y]
tform3 = fitgeotform2d(movingPoints, fixedPoints, 'similarity'); %once

refDMD = imref2d(mrsize_xy);
refDMDbig = imref2d(usFactor * mrsize_xy);
refDMDbig.XWorldLimits = refDMD.XWorldLimits;
refDMDbig.YWorldLimits = refDMD.YWorldLimits;

%% texture mapping from 2D to 3D in Allen CCF (Kim) space
Vus = imresize3(V, usFactor, 'linear');%once
[surfDepth] = getSurfaceData3(Vus, 'last', 0);%once
sizeVus = size(Vus);

% testimg = ones(imagesize);
% testimg_reg = imwarp(testimg,tform3,'linear', ...
%     'OutputView', refDMDbig, 'FillValues',0); %all
% TexVol = paintSurfaceToVolume(surfDepth, testimg_reg, sizeVus);%all
