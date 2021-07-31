addpath('C:\Users\Sato Lab\Documents\MATLAB\Suite2Ptool');


load('G:\Data5\MO137\Imaging\20190522\TurboRegImage\MO137_20190522_TurboReg__neuropil.ext', '-mat')

load('K:\Suite2Pdata\MO137\20190522\F_MO137_20190522_plane1_proc.mat',...
    'dat');
ImageInfos_s2p = Marius2Takashi(dat);

%% check1: ROI
figure;
for jj = 1:2   
    switch jj
        case 1
            thisImageInfos = ImageInfos; %MUSC
        case 2
            thisImageInfos = ImageInfos_s2p; %suite2p 
    end
    
    subplot(1,2,jj);
    nCell = size(thisImageInfos.filterMatrix,1);
    imagesc(dat.ops.mimg);
    switch jj
        case 1
            title('MUSC');
        case 2
            title('Suite2p');
    end
    
    hold on;
    for iCell = 1:nCell
        contour(squeeze(thisImageInfos.filterMatrix(iCell,:,:)),'color','w');
%         pause(0.1);
%         contour(squeeze(thisImageInfos.filterMatrix_neuropil(iCell,:,:)),'color','g');
%         plot(thisImageInfos.centerPosition(iCell,1), thisImageInfos.centerPosition(iCell,2), 'rx');
    end
    axis equal tight
end

saveas(gcf,'ROI.fig');


%% check2: trace
%location of target cell to show trace
xc=360;
yc=120;

xc=290;
yc=420;

figure;
for jj= 1:2
    switch jj
        case 1
            thisImageInfos = ImageInfos; %MUSC
        case 2
            thisImageInfos = ImageInfos_s2p; %suite2p
    end
    
    [~,thisCell] = min((thisImageInfos.centerPosition(:,1)-xc).^2 ...
        + (thisImageInfos.centerPosition(:,2)-yc).^2);
    
    for iTr=1:7
        thisImageData = thisImageInfos.ImageData(iTr);
        nFrames = size(thisImageData.Fvalue(thisCell,:),2);
        taxis = thisImageData.relativeTime(thisCell,1:nFrames);
       
        subplot(7,2,2*(iTr-1)+jj);
        plot(taxis, thisImageData.Fvalue(thisCell,:));
        hold on
        plot(taxis, thisImageData.Fvalue_neuropil(thisCell,:));
        title(['trialNo: ' num2str(thisImageData.trialNo) ...
            ', sliceNoList ' num2str(thisImageData.sliceNoList(thisCell)) ...
            ', validityList: ' num2str(thisImageData.validityList(thisCell))])
        grid minor
    end
    
end
linkaxes;
marginplots;

saveas(gcf,'traces.fig');

