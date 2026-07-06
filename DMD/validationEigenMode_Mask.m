

resultServer = '/mnt/dshi0006_market';

modeIdx = [2 3 4 20 50 100];
nsites = 1;
for ii =1:numel(modeIdx)
    imgPrefix_mode = ['whole_mode=' num2str(modeIdx(ii))];
    imgPrefix_mask = ['mask_whole_mode=' num2str(modeIdx(ii)) '_numsites=' num2str(nsites)];

    mode = load(fullfile(resultServer, 'DMD images',imgPrefix_mode, [imgPrefix_mode '_Confucious.mat']));
    mask = load(fullfile(resultServer, 'DMD images',imgPrefix_mask, [imgPrefix_mask '_Confucious.mat']));

    ax(1)=subplot(121);
    imagesc(mode.image4OI(:,:,1));axis equal tight;title('positive');hold on;
    contour(mask.image4OI(:,:,1),1,'r');
    ax(2)=subplot(122);
    imagesc(mode.image4OI(:,:,2));axis equal tight;title('negative');hold on;
    contour(mask.image4OI(:,:,2),1,'r');
    linkcaxes(ax);
    mcolorbar;

    screen2png(['validateMode_Mask_' imgPrefix_mode]);
    close all
end