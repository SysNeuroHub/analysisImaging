function refImage = getReference(subjectName, tgt)
% return image of reference image of a subject
%if tgt = 'tube', only show reference tubes
%if tgt = 'brain', show both brain and tubes

%UNDER CONSTRUCTION

%load widefield image under amber LED exposure
resultServer = '/mnt/dshi0006_market';
load(fullfile(resultServer, 'Subjects', subjectName, 'MR2DMDresult','Atlas_reg_info.mat'),  ...
    'image2','image2th');

%% camera image info
refDir = '/home/daisuke/Documents/git/analysisImaging/DMD/references';
refdate = '20260214';
load(fullfile(refDir, ['camImg_' refdate]),'camImg');
DMDprojectionZone = getDMDprojectionZone(camImg);

if strcmp(tgt, 'tube')
    thisImage = normalize_prctile(image2,[image2th,100]);
elseif strcmp(tgt, 'brain')
    thisImage = log(image2);
end

f = figure('visible','off');
imagesc(thisImage); colormap(gray); axis equal tight; hold on
contour(DMDprojectionZone,'r');
plot(camImg.bregmapix(2), camImg.bregmapix(1), 'kx')
% plot(camImg.lambdapix(2), camImg.lambdapix(1), 'yx'); %DO NOT USE
%screen2png(fullfile(imgDir, [imgPrefix '_' subject '_all_wCCF']));
exportPng4DMD('refimg_tmp', f, 0);
refImage = imread('refimg_tmp.png');
delete('refImg_tmp.png');

