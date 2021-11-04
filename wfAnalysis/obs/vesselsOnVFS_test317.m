%% draw signmap and blood vessels
%TODO
% use a picture focused on brain surface
% independent control of crange between blue and purple channels
% add scale bar and 1mm grid
% convert pref phase to visual angle

expt.subject = 'TIGRE2GCaMP6s_317';
expt.expDate = '2020-10-15_1';
expt.expNum = 1;

resultSaveDir = fullfile('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\AnalysisResult\wf',...
    [expt.subject '_' expt.expDate '_' num2str(expt.expNum)]);

load(fullfile(resultSaveDir,[figname 'signMap']), 'signMap_median','xMapm','yMapm','resizeS');
signMap_median = imresize(signMap_median, 1/resizeS);
xMapm = imresize(xMapm, 1/resizeS);
yMapm = imresize(yMapm, 1/resizeS);

movieSuffix = {'blue','purple'};
crange = [1 80];
mimg_norm = [];
for isuffix = length(movieSuffix)
    params.movieSuffix = movieSuffix{isuffix};
   [~, ~, ~, mimg, mask] = quickLoadUVt(expPath, nSV, saveVpath, params);
   
    %normalize 0-1 for transparency
    range_p = prctile(mimg(mask),crange);
    mimg_n = 1/diff(range_p)*(mimg - range_p(1));
    mimg_n(mimg_n<0)=0;
    mimg_n(mimg_n>1)=1;
    
    mimg_norm(:,:,isuffix) = mimg_n;

    
    
    %% visuzlize
    figure('position',[0 0 1500 400]);
    whitebg(gcf,'k'); %so that the vessels look darl
    
    ax(1)=subplot(131);imagesc(xMapm,'alphadata',mimg_norm(:,:,isuffix));axis equal tight;
    caxis(prctile(xMapm(:),[1 99]));mcolorbar(gca,.5);title('pref x');
    axis off
    
    ax(2)=subplot(132);imagesc(yMapm,'alphadata',mimg_norm(:,:,isuffix));axis equal tight;
    caxis(prctile(yMapm(:),[1 99]));mcolorbar(gca,.5);title('pref y');
    axis off
    
    ax(3)=subplot(133);imagesc(signMap_median,'alphadata',mimg_norm(:,:,isuffix));axis equal tight;
    caxis([-1 1]);mcolorbar(gca,.5);title('field sign');
    axis off
    colormap(cool);
    
    linkaxes(ax,'xy');
    set(gcf,'color','k');
    
    screen2png(fullfile(resultSaveDir,[figname '_vessel' movieSuffix{isuffix}]), gcf);
end
