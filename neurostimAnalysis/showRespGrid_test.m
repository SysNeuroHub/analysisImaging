%parameter for analysis
winResp = [2 5];%[3 6];
winBase = [winSamps(1) 0];


%input ... result of eventLockedAvg
%winSamps
%avgPeriEventStack
eLabels = unique(stimLabels);


[~,framesResp(1)] = min(abs(winSamps-winResp(1)));
[~,framesResp(2)] = min(abs(winSamps-winResp(2)));
[~,framesBase(1)] = min(abs(winSamps-winBase(1)));
[~,framesBase(2)] = min(abs(winSamps-winBase(2)));

mEvk = mean(avgPeriEventStack(:,:,framesResp(1):framesResp(2),:),3);
mBase = mean(avgPeriEventStack(:,:,framesBase(1):framesBase(2),:),3);

mResp = squeeze(1e2*(mEvk - mBase)./mBase); %dF/F [%] 

xpos = real(eLabels);
ypos = imag(eLabels);

xparam = unique(xpos);
yparam = unique(ypos);

nXgrid = length(xparam);
nYgrid = length(yparam);

ax=[];
for ii = 1:nXgrid*nYgrid
    fp = find(xparam == xpos(ii)) + (find(yparam == ypos(ii))-1)*nXgrid; %figure position
    thisCond = find(eLabels == xpos(ii) + 1j*ypos(ii));
    ax(ii) = subplot(nYgrid,nXgrid,fp);
    imagesc(mResp(:,:,thisCond));
    axis equal tight;
    title(['x: ' num2str(xpos(ii)) ', y: ' num2str(ypos(ii))]);
    grid minor
    caxis(prctile(mResp(:),[0 99]));
end
linkaxes(ax(:));
subplot(nYgrid,nXgrid,nXgrid*nYgrid);
mcolorbar;
