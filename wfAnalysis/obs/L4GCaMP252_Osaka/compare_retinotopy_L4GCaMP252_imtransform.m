load('obs_to_new','t_concord');


load('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\AnalysisResult\wf\L4GCaMP6s_252_2021-01-09_2_1\2021-01-09_2_1_L4GCaMP6s_252_stimKalatsky_DS_blue_correctedsignMap.mat')

 blue = loadTiffStack('\\vault-v2.erc.monash.edu.au\MNHS-dshi0006\Subjects\L4GCaMP6s_252\2021-01-09\blue.tif',...
     'imread',0);
 blue = imresize(double(blue),0.5);
 imsize_new = size(blue);
 
 %common between the two exps
 roix = 140:215;
 roiy = 25:110;

 figure(1);
h1=subplot(131);
mimg_r=log(double(flipud(blue(roiy,roix))'));
imagesc(mimg_r);colormap(h1,gray);
caxis(prctile(mimg_r(:),[10 99]));
axis equal tight off;

h2=subplot(132);
thisimage=flipud(signMap_median(roiy,roix))';
imagesc(thisimage);colormap(h2,jet);
caxis([-1 1]);
axis equal tight off;
mcolorbar(gca,.5);

h3=subplot(133);
thisimage=flipud(xMapm(roiy,roix))';
imagesc(thisimage);colormap(h3,jet);
caxis(prctile(thisimage(:),[10 90]));
axis equal tight off;
screen2png('new');


%% old
load('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\AnalysisResult\wf\L4GCaMP6s_252_2020-11-05_1_1\2020-11-05_1_1_L4GCaMP6s_252_stimKalatsky_DS_blue_correctedsignMap.mat');
blue = loadTiffStack('\\vault-v2.erc.monash.edu.au\MNHS-dshi0006\Subjects\L4GCaMP6s_252\2020-11-05\blue focus.tif',...
    'imread',0);
blue = imresize(double(blue),0.5);

%align to the new one
blue = imtransform(blue, t_concord,'XData',[1 imsize_new(2)],'YData',[1 imsize_new(1)]);
signMap_median = imtransform(signMap_median, t_concord,'XData',[1 imsize_new(2)],'YData',[1 imsize_new(1)]);
xMapm = imtransform(xMapm,t_concord,'XData',[1 imsize_new(2)],'YData',[1 imsize_new(1)]);

%  roix = 100:(100+75);
%  roiy = 30:(30+85);
 
figure(2);
h1=subplot(131);
mimg_r=log(double(flipud(blue(roiy,roix))'));
imagesc(mimg_r);colormap(h1,gray);
caxis(prctile(mimg_r(:),[10 99]));
axis equal tight off;

h2=subplot(132);
thisimage=flipud(signMap_median(roiy,roix))';
imagesc(thisimage);colormap(h2,jet);
%caxis(prctile(thisimage(:),[10 90]));
caxis([-1 1]);
axis equal tight off;
mcolorbar(gca,.5);
h3=subplot(133);
thisimage=flipud(xMapm(roiy,roix))';
imagesc(thisimage);colormap(h3,jet);
caxis(prctile(thisimage(:),[10 90]));
axis equal tight off;
screen2png('obs');
