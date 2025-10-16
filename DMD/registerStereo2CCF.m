function [tform3, surfDepth, Vusinfo] = registerStereo2CCF(bregma, MmPerPixel, ...
    path2Brain_template, usFactor)
% [tform3, surfDepth, Vusinfo] = registerStereo2CCF(bregma, MmPerPixel, ...
%     path2Brain_template, usFactor)
% computes 2D transformation from stereotaxi coordinate to CCF looked from
% above, with upsampling factor usFactor

if isempty(usFactor)
    usFactor = 5;
end

MmPerPixel_mr = 0.1; %[mm]


%V = flip(niftiread(path2Brain_template)>0, 3);
V = afterniftiread(niftiread(path2Brain_template)>0);
V_info = niftiinfo(path2Brain_template);
mrsize_xy = [size(V,3), size(V,1)]; %[x y]


bregmak_c = allenCCFbregma()*MmPerPixel_mr;
bregmak = [bregmak_c(1) bregmak_c(3)]+0.5; %[y x]
movingPoints = [bregma(2) bregma(1); bregma(2)+1/MmPerPixel bregma(1)]; %[x y]
fixedPoints = [bregmak(2) bregmak(1); bregmak(2)+1/MmPerPixel_mr bregmak(1)]; %[x y]
tform3 = fitgeotform2d(movingPoints, fixedPoints, 'similarity'); %once

refDMD = imref2d(mrsize_xy);
refDMDbig = imref2d(usFactor * mrsize_xy);
refDMDbig.XWorldLimits = refDMD.XWorldLimits;
refDMDbig.YWorldLimits = refDMD.YWorldLimits;


% Vus = imresize3(V, usFactor, 'linear');
% [surfDepth] = getSurfaceData(Vus, 'last', 0.5);

V = smooth3(V, 'gaussian', [7 7 7], 3); %important for projection mapping
Vus = imresize3(V, usFactor, 'linear');
[surfDepth] = getSurfaceData(Vus, 'last', 0.5);
%sizeVus = size(Vus);

Vusinfo = V_info;
Vusinfo.Datatype = 'double';
Vusinfo.Transform.T(1:3,1:3) = Vusinfo.Transform.T(1:3,1:3)/usFactor;
Vusinfo.PixelDimensions = Vusinfo.PixelDimensions/usFactor;
Vusinfo.ImageSize = size(Vus);


% testimg = ones(imagesize);
% testimg_reg = imwarp(testimg,tform3,'linear', ...
%     'OutputView', refDMDbig, 'FillValues',0); %all
% TexVol = paintSurfaceToVolume(surfDepth, testimg_reg, sizeVus);%all
