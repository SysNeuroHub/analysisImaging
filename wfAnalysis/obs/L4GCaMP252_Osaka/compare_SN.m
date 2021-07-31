load('obs_to_new','t_concord');


for dd=1:2
    switch dd
        case 1 %new
            load('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\AnalysisResult\wf\L4GCaMP6s_252_2021-01-09_2_1\2021-01-09_2_1_L4GCaMP6s_252_stimKalatsky_DS_blue_correctedsignMap.mat')
            gavgABSMAPS=squeeze(mean(mean(ABSMAPS,3),4));
            [a,ind]=max(gavgABSMAPS(:));
            [peak_y, peak_x] = ind2sub(size(gavgABSMAPS),ind);
            
            imsize_new = size(gavgABSMAPS);
        case 2 %old
            load('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\AnalysisResult\wf\L4GCaMP6s_252_2020-11-05_1_1\2020-11-05_1_1_L4GCaMP6s_252_stimKalatsky_DS_blue_correctedsignMap.mat');
            ABSMAPS = imtransform(ABSMAPS, t_concord, 'XData',[1 imsize_new(2)],'YData',[1 imsize_new(1)]);
            ABSMAPS_u = imtransform(ABSMAPS_u, t_concord, 'XData',[1 imsize_new(2)],'YData',[1 imsize_new(1)]);
            ABSMAPS_d = imtransform(ABSMAPS_d, t_concord, 'XData',[1 imsize_new(2)],'YData',[1 imsize_new(1)]);
            
    end
%     S = ABSMAPS(peak_y,peak_x,1:4,:);
%     N = 0.5*(ABSMAPS_u(peak_y,peak_x,1:4,:) + ABSMAPS_d(peak_y,peak_x,1:4,:));
%     meanS= mean(S(:));
%     meanN = mean(N(:));
    S = squeeze(mean(mean(ABSMAPS(:,:,1:4,:),3),4));
    N = squeeze(mean(mean(0.5*(ABSMAPS_u(:,:,1:4,:) + ABSMAPS_d(:,:,1:4,:)),3),4));

    SN(:,:,dd) = S./N;
  
end

for dd = 1:2
      subplot(2,1,dd);
    imagesc(squeeze(SN(:,:,dd)));
    caxis(prctile(SN(:),[5 95]));
    axis equal tight
end

 roix = 140:215;
 roiy = 25:110;

 sn_part = SN(roiy,roix,:);
 sn_part = reshape(sn_part,length(roix)*length(roiy),2);
 sn_part = fliplr(sn_part);
 
%  for dd=1:2
%  scatter(dd*ones(length(sn_part),1),sn_part(:,dd),[],...
%      ([.7 .7 .7].*ones(length(sn_part),1)),'.');
%  hold on;
%  end
 boxplot(sn_part);
 set(gca,'xticklabel',{'week 1', 'week 10'});
 screen2png('compare_SN');
 