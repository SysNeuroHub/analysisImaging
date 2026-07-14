function f = showSingleTracesSVD(U,V,t,xx,yy,icond,expt,stimSequence, respWin)
%f = showSingleTracesSVD(U,V,t,expt,stimSequence, respWin)

% respWin = [-0.5 1.5];
[avgPeriEventV, winSamps, periEventV, sortLabels] = ...
    eventLockedAvg(V, t, expt.stimTimes.onset, stimSequence.seq, respWin);
preIdx = find(winSamps<0);
f = figure('position',[0 0 900 420]);
%yy=123; xx = 69;
%icond = 17;
these = find(sortLabels == icond);
trace = [];
for ii = 1:numel(these)
    trace(:,ii)=svdFrameReconstruct(U(yy,xx,:),squeeze(periEventV(these(ii),:,:)));
end

subplot(121);
imagesc(winSamps, 1:numel(these), trace');vline(0);
climF = round(prctile(trace(:),[1 99]));
title(['F, y=' num2str(yy) ', x=' num2str(xx) ', icond:' num2str(icond)]);
caxis(climF);
mcolorbar(gca, .5);
xlabel('Time from laser onset [s]');
ylabel('trial');


subplot(122);
FF0 = trace'-mean(trace(preIdx,:))';
climFF0 = round(prctile(abs(FF0(:)),99));
imagesc(winSamps, 1:numel(these), FF0);vline(0);
title('F-F0')
caxis([-climFF0 climFF0])
mcolorbar(gca, .5);