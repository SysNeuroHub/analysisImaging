%% save data
% setPath_analysisImaging;
% 
% expt.subject = 'Tiberius';
% expt.expDate = '2026-07-10_1';
% expt.expNum = 1;
% params.polyScanPro = false;
% params.bklightCtrl = 0;
% params.movieSuffix = 'amber';
% params.useCorrected = 1;
% params.nSV = 1000;
% params.resizeS = 0.5; %spatial rescaling factor
%
% packed = packngo(expt, params);
% save('Tiberius_2026-07-10_1','packed','params');


%% load data
load('Tiberius_2026-07-10_1','packed','params');
expt = packed.expt;


%% event-triggered analysis
respWin = [-0.2 0.6]; %[s]
istim = 11; 

% brain signals triggered by DMD stimulus onset
[avgPeriEventV, winSamps, periEventV, sortLabels] = ...
    eventLockedAvg(packed.V, packed.t, expt.laserTimes.onset, expt.stimSequence.seq, respWin);
avgPeriEventV = avgPeriEventV - mean(avgPeriEventV(:,:,winSamps<0),3); % subtract pre-DMD-stimulus period

% DMD images triggered by DMD stimulus onset ... NEEDS IMPROVEMENT
DMDimgSequence = getDMD4OI_eventLockedAvg(packed.DMDimg, ...
    expt.laserAmp, expt.DMDidx, packed.t, expt.stimSequence, ...
    expt.laserTimes, respWin);


%% response to 11th DMD stimuli, averaged across repeats

avgResponse = svdFrameReconstruct(packed.U, squeeze(avgPeriEventV(istim,:,:)));

ax_DMD = axes;
ax_brain = axes;
nPanels = floor(.5*numel(winSamps));
for ii = 1:nPanels
    tidx = 2*ii;
    ax_DMD(ii) = subplot(2,nPanels, ii);
    imagesc(squeeze(DMDimgSequence(:,:,tidx,istim)),'AlphaData',1-packed.ccf); axis square equal off;
    title(['t=' num2str(winSamps(tidx)) '[s]']);
    
    ax_brain(ii) = subplot(2, nPanels, ii+nPanels);
    imagesc(squeeze(avgResponse(:,:,tidx)),'AlphaData',1-packed.ccf);  axis square equal off;
end

colorbar(ax_DMD(end));
colorbar(ax_brain(end));
linkcaxes(ax_DMD);
linkcaxes(ax_brain);


%% cf. GUI to show event-triggered trace of all DMD stimuli
% pixelTuningCurveViewerSVD(packed.U, packed.V, packed.t, expt.laserTimes.onset, ...
%     expt.stimSequence.seq, respWin, 1);


