function [TexImgwarpedtoDMD, TexImgwarped] = applyCCF2DMD(VCCF, tform_T2OI, tform_OIDMD, OIsize, MmPerPixel_oi, MROIDMDsubjectDir, autoTform, extractSurface)
% TexImgwarpedtoDMD = applyStereo2DMD(VCCF, mrangle, tform_T2OI, tform_OIDMD, OIsize, MmPerPixel_oi, MROIDMDsubjectDir)
% transforms images defined in CCF volume to DMD space
% input must be [0 1]
%
% need
% TexVol_T2.nii (in Stereo2T2.m)
% TexVol.nii (in Stereo2T2.m)
% Brain_template.nii
% T2w_brain.nii
% 
% created from applyStereo2DMD
if nargin < 8
    extractSurface = 1;
end

if nargin < 7
    autoTform = 1;
end
if nargin < 5
    MmPerPixel_oi = [];
end

verbose = 0;

%OIsize = [1080 1080];  %[y x]
%MmPerPixel_oi = 0.0104;
DMDsize = [500 800]; %[y x]
MmPerPixel_mr = 0.1; %[mm]
% usFactor = 5;
path2Brain_template =  fullfile(MROIDMDsubjectDir,'Brain_template.nii');

% tform_T2OI_us = tform_T2OI;
% if isfield(tform_T2OI, 'Scale')
%     tform_T2OI_us.Scale = tform_T2OI.Scale/usFactor;
% end

% % if isempty(MmPerPixel_oi)
% %     MmPerPixel_oi = MmPerPixel_img;
% % end
% % 
% % if islogical(images)
% %     error('Input image must be numeric, not logical');
% % end
% % 
% % if max(images(:)>1)||min(images(:)<0)
% %     error('Input image must be [0 1]');
% % end
nVols = size(VCCF,4);

% % disp('Running registerStereo2CCF...');
% % [tform3, surfDepth, Vusinfo] = registerStereo2CCF(bregma, MmPerPixel_img, ...
% %     path2Brain_template, usFactor);
% % 
% % %for AtlasTexVol_to_T2.sh
% % niftiwrite_us(fullfile(MROIDMDsubjectDir,'T2w_brain.nii'), usFactor);

% extracted and modified from register2Stereo2CCF
V = afterniftiread(niftiread(path2Brain_template)>0);
V = smooth3(V, 'gaussian', [7 7 7], 3); 

if extractSurface
    [xx,yy,zz] = ndgrid(-2:2);
    nhood = sqrt(xx.^2 + yy.^2 + zz.^2) <= 2.0;
    Verode = imerode(V, nhood); %important to get surfDepth
    [surfDepth] = getSurfaceData(Verode, 'last', 0.95);
else
    surfDepth = [];
end
V_info = niftiinfo(path2Brain_template);
V_info.Datatype = 'double';


currentDir = pwd;
cd(MROIDMDsubjectDir);

if extractSurface
    surfaceT2 = CCF2T2(ones(size(VCCF(:,:,:,1))), V_info, surfDepth);
end

TexImgwarped = zeros(OIsize(1), OIsize(2), nVols);
TexImgwarpedtoDMD = zeros(DMDsize(1), DMDsize(2), nVols);
for ii = 1:nVols
    disp(['applyCCF2DMD: ' num2str(ii) '/' num2str(nVols)]);

    
    TexVol_T2 = CCF2T2(VCCF(:,:,:,ii), V_info, surfDepth);

    if extractSurface
        TexImg = squeeze(nanmean(TexVol_T2./surfaceT2,3));
    else
        TexImg = squeeze(nanmean(TexVol_T2,3));
    end

       if verbose
        subplot(1,3,1); imagesc(TexImg);axis equal tight; 
        title('CCF warped to upscaled T2');
    end

    %% project to OI
    if autoTform
        fixedRef  = imref2d(OIsize, MmPerPixel_oi, MmPerPixel_oi);  % example pixel sizes in mm
        movingRef = imref2d(size(TexImg), MmPerPixel_mr,  MmPerPixel_mr);
        TexImgwarped(:,:,ii) = imwarp(TexImg,movingRef,tform_T2OI,'linear','OutputView',fixedRef);
    else
        TexImgwarped(:,:,ii) = imwarp(TexImg,tform_T2OI,'linear','OutputView',imref2d(OIsize));
    end

    if verbose
        subplot(1,3,2); imagesc(TexImgwarped(:,:,ii) );axis equal tight; 
        title('T2 warped to OI'); 
    end

    %% project to DMD
    TexImgwarpedtoDMD(:,:,ii) = imwarp(TexImgwarped(:,:,ii) ,tform_OIDMD,'linear', ...
        'OutputView',imref2d(DMDsize));
    if verbose
        subplot(1,3,3); imagesc(TexImgwarpedtoDMD);axis equal tight; 
        title('OI warped to DMD');

        fig = gcf;
        %caxes(fig, prctile(images(:), [0 100]));
    end
end
delete('TexVol.nii', 'TexVol_T2.nii');
cd(currentDir);

