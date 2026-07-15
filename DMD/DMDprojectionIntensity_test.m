%uiopen('/mnt/dshi0006_market/tmp/allmirroron20260708.tif',1);
% thisImage = allmirroron20260708;

%from getDMDprojectionZone
% imgDir = '/home/daisuke/Documents/git/analysisImaging/DMD/references/DMDprojectionZone';
% latestImageFile = 'DMDprojectionZone20260214.tif';
% DMDfullOn = imread(fullfile(imgDir, latestImageFile));
% thisImage = resizeCropFullImg(DMDfullOn, camImg);
refdate = '20260214';
refDir = '/home/daisuke/Documents/git/analysisImaging/DMD/references';
load(fullfile(refDir, ['camImg_' refdate]),'camImg');
[DMDprojectionZone, thisImage] = getDMDprojectionZone(camImg);

%% raw image
[nRows, nCols] = size(thisImage);
subplot(2,2,1);
imagesc(thisImage);mcolorbar
vline(nCols/2);hline(nRows/2);
subplot(222);
plot(thisImage(:,nCols/2),1:nRows);axis tight ij;
subplot(223);
plot(thisImage(nRows/2,:));axis tight;

grad = imgradient(thisImage);

proj = (grad.*DMDprojectionZone);
ngInd = find(abs(proj - mean(proj(:))) > 3*std(proj(:)));
thisImage_i = double(thisImage);%.*DMDprojectionZone;
thisImage_i(ngInd) = nan;
thisImage_i = fillmissing2(thisImage_i, "movmedian",10);
thisImage_i(~DMDprojectionZone) = nan;

thisImage_is = smoothdata(thisImage_i);
%thisImage_is = imgaussfilt(thisImage_i,'FilterSize',9); 

figure
ax(1)=subplot(131);
imagesc(thisImage);axis equal tight;mcolorbar;title('original');
ax(2)=subplot(132);
imagesc(thisImage_i);axis equal tight;mcolorbar;title('interpolated');
ax(3)=subplot(133);
imagesc(thisImage_is);axis equal tight;mcolorbar;title('smoothed');
linkcaxes(ax);
