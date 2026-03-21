function TexVol_T2 = CCF2T2(VCCF, V_info, surfDepth)
% TexVol_T2 = Stereo2T2(VCCF, V_info, surfDepth)
% projectionmaps images defined in CCF volume then to T2 space
% need
% TexVol_T2.nii
% TexVol.nii

% created from Stereo2T2.m


verbose = 0;
% sizeVus = Vusinfo.ImageSize;

% if verbose
%     fig = figure;
%     subplot(1,3,1); imagesc(imageStereo); hold on; axis equal tight;
%     title('image in stereotaxic coords');
% end
% 
% %% warp from stereotaxic to CCF (2D)
% refDMD = imref2d(round([sizeVus(1), sizeVus(2)]/usFactor));
% refDMDbig = imref2d([sizeVus(1), sizeVus(2)]);
% refDMDbig.XWorldLimits = refDMD.XWorldLimits;
% refDMDbig.YWorldLimits = refDMD.YWorldLimits;

% img_reg = imwarp(imageStereo,tform3,'linear', ...
%     'OutputView', refDMDbig, 'FillValues',0);
% 
% if verbose
%     subplot(1,3,2); imagesc(img_reg); axis equal tight;
%     title('stereo image registered to upscaled CCF');
% end

%% extract pixel values on brain surface
if nargin > 2 && ~isempty(surfDepth)
    TexVol_CCF = beforeniftiwrite(extractSurfaceValues(surfDepth, VCCF));
else
    TexVol_CCF = beforeniftiwrite(VCCF);
end

if verbose
    subplot(1,3,3); isosurface(TexVol_CCF); axis equal tight;
    f = gcf;
    f.CurrentAxes.YDir = 'Reverse';
    ax = gca; set(gca, 'view', [120 20]);
    title('Projection mapped to upscaled CCF');
end

% V_info4saving = V_info;
% V_info4saving.ImageSize = [V_info.ImageSize(2) V_info.ImageSize(3) V_info.ImageSize(1)];
niftiwrite(TexVol_CCF, 'TexVol.nii', V_info);

%% transform from CCF to individual T2
cmdStr = [fullfile(fileparts(mfilename('fullpath')), 'AtlasTexVol_to_T2.sh') ' ' 'TexVol.nii' ' ' 'T2w_brain.nii' ' ' pwd];
system(cmdStr);
%expected to exist under pwd: Tem_to_T21Warp.nii.gz, Tem_to_T20GenericAffine.mat
%output: TexVol_T2.nii

TexVol_T2 = afterniftiread(niftiread('TexVol_T2.nii'));
