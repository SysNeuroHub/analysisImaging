function [Ub, Vcorr, tb, mimgB] = quickHemoCorrect(expPath, savePath, nSV, ...
    hemoFreq, pixSpace, suffix_fluo, suffix_hemo)
% [Ub, Vcorr, tb, mimgB] = quickHemoCorrect(expPath, savePath, nSV, hemoFreq, pixSpace)
% loads SVD of suffix_fluo(blue) and suffix_hemo(purple),
% transforms V of blue so that timestamp is common between blue and purple (alignTimeStamp.m),
% transforms V of purple so that U is common between blue and purple (ChangeU.m),
% returns V of blue after hemodynamic correction (HemoCorrectLocal.m)
%
% Inputs:
%   expPath: path to where original SVD is saved
%   savePath: path to where created SVD will be saved
%   nSV: number of SVs to compute hemodynamics correction
%   hemoFreq: frequency range at which the correction factor is computed [low high] (Hz)
%   pixSpace: subsampling factor (in pixel) to compute the correction factor
%
% Outputs:
%   Ub: U of blue channel (identical to before hemodynamic correction)
%   Vcorr: V of blue channel after hemodynamic correction
%   tb: time stamp of blue channel (identical to before hemodynamic correction)
%   mimgB: mean image of blue channel (identical to before hemodynamic correction)

detrendFilt = 1; %11/6/20 if true apply detrend&low cut fltering to Vb & Vp

if nargin < 4
    hemoFreq = [10 13]; % frequency to look for hemo signals, can also be heartbeat [9,13]
end
if nargin < 5
    pixSpace = 3; % Something about subsampling the image to correct hemo
end
if nargin < 6
    suffix_fluo = 'amber';
end
if nargin < 7
    suffix_hemo = 'red';
end

expRoot = fileparts(expPath);
fprintf(1, ['loading ' suffix_fluo '\n'])
try
    Ub = readUfromNPY(fullfile(expRoot, ['svdSpatialComponents_' suffix_fluo '.npy']), nSV);
    mimgB = readNPY(fullfile(expRoot, ['meanImage_' suffix_fluo '.npy' ]));
    load(fullfile(expRoot, ['dataSummary_' suffix_fluo '.mat']));
    DSb = dataSummary;
catch err
    Ub = readUfromNPY(fullfile(expPath, ['svdSpatialComponents_' suffix_fluo '.npy']), nSV);
    mimgB = readNPY(fullfile(expPath, ['meanImage_' suffix_fluo '.npy']));
    load(fullfile(expPath, ['dataSummary_' suffix_fluo '.mat']));
    DSb = dataSummary;
end
Vb = readVfromNPY(fullfile(expPath, ['svdTemporalComponents_' suffix_fluo '.npy']), nSV);
tb = readNPY(fullfile(expPath, ['svdTemporalComponents_' suffix_fluo '.timestamps.npy']));

fprintf(1, ['loading ' suffix_hemo '\n'])
try
    Up = readUfromNPY(fullfile(expRoot, ['svdSpatialComponents_' suffix_hemo '.npy']), nSV);
    mimgP = readNPY(fullfile(expRoot, ['meanImage_' suffix_hemo '.npy']));
    load(fullfile(expRoot, ['dataSummary_' suffix_hemo '.mat']));
    DSp = dataSummary;
catch err
    Up = readUfromNPY(fullfile(expPath, ['svdSpatialComponents_' suffix_hemo '.npy']), nSV);
    mimgP = readNPY(fullfile(expPath, ['meanImage_' suffix_hemo '.npy']));
    load(fullfile(expPath, ['dataSummary_' suffix_hemo '.mat']));
    DSp = dataSummary;
end
Vp = readVfromNPY(fullfile(expPath, ['svdTemporalComponents_' suffix_hemo '.npy']), nSV);
tp = readNPY(fullfile(expPath, ['svdTemporalComponents_' suffix_hemo '.timestamps.npy']));

if size(Vb,2)>size(Vp,2)
    % can be an extra blue frame - need same number
    Vb = Vb(:,1:end-1);
    tb = tb(1:end-1);%15/6/20
end

Fs = 1/mean(diff(tb));

[nRows, nCols, nFrames] = size(Ub); 
yidx = round(nRows/4):round(nRows*3/4);%701:725;
xidx = round(nCols/4):round(nCols*3/4);%301:325;
roi = zeros(size(Ub));
roi(yidx,xidx,:) = 1;

if strcmp(hemoFreq,'auto')
    ROItrace = squeeze(svdFrameReconstruct(mean(mean(Up.*roi)),Vp));

    L = length(ROItrace);
    NFFT = 2^nextpow2(L);
    [Pxx,F] = pwelch(ROItrace-mode(ROItrace),[],[],NFFT,Fs);
    %    plot(F, 10*log10(Pxx));

    % rmpath(genpath('/home/daisuke/Documents/git/encoding'));
    % f = fit(F,10*log10(Pxx),'exp1');
    % plot(f,F,10*log10(Pxx));
    theseFreqInds = find(F>1);
    [~, maxInds] = max(Pxx(theseFreqInds));
    hemoFreq = [max(1,F(theseFreqInds(maxInds))-2) min(F(theseFreqInds(maxInds))+2, 0.99*Fs/2)];%15/10/20
    disp(['Hemofreq detected as: ' num2str(hemoFreq)]);
end


%% detrend amber and red in original space
[imageMeans_pc,base_p,param_p,diag_p] = bleachCorrectionRobustNorm(DSb.imageMeans);
corrected_p = bleachCorrectionRobustNorm(svdFrameReconstruct(Up,Vp), param_p.beta); %too heavy

%% re-SVD of amber and red
    [Ub, Svb, Vb, totalVarb] = get_svdcomps(svdOps);
    [Up, Svp, Vp, totalVarp] = get_svdcomps(svdOps);

% if detrendFilt
% 
%      Vb = detrend(Vb', 'linear')'+mean(Vb, 2);
%      Vp = detrend(Vp', 'linear')'+mean(Vp, 2);
% end

% align blue/purple in time, apply correction
Vbs = alignTimeStamp(Vb, tb, tp, 'interp'); %10/6/20 as in mergeST.m


VpNewU = ChangeU(Up, Vp, Ub); % puts hemo-color V into signal-color U space
[Vcorr,~,ScaleFactor] = HemoCorrectLocal(Ub, Vbs, VpNewU, Fs, hemoFreq, pixSpace); % correct for hemodynamics

% svdViewer(Ub, DSb.Sv, Vb, Fs) %blue before correction fig1
% svdViewer(Up, DSp.Sv, Vp, Fs) %purple fig2
% svdViewer(Ub, DSb.Sv, Vcorr, Fs) %blue after correction fig3
% newfig = mergefigs([1 2 3]);%merge figures produced in quickHemoCorrect
% set(newfig,'position',[1 1 1920 1200]);
% screen2png(fullfile(saveVpath,'HemoCorrect_sanityCheckFigs'));


ROItrace_b = squeeze(svdFrameReconstruct(nanmean(nanmean(Ub.*roi)), Vb)); %squeeze(nanmean(nanmean(svdFrameReconstruct(Ub,Vb))));
ROItrace_p = squeeze(svdFrameReconstruct(nanmean(nanmean(Up.*roi)), Vp));
ROItrace_corr = squeeze(svdFrameReconstruct(nanmean(nanmean(Ub.*roi)), Vcorr));


fig=figure('position',[0 0 1900 1200]);
ax(1)=subplot(411);
plot(tb, ROItrace_b, 'Color',[1 .5 0]); hold on;
plot(tp, ROItrace_corr,'b');
grid minor
legend(suffix_fluo, [suffix_fluo ' corrected']);

ax(2)=subplot(412);
plot(tp, ROItrace_p, 'r');
grid minor
legend(suffix_hemo);
linkaxes(ax,'x');
xlim([10 40]);

subplot(223);
[Pxxb, Fb] = myTimePowerSpectrum(ROItrace_b, Fs);
[Pxxp, Fp] = myTimePowerSpectrum(ROItrace_p, Fs);
[Pxxcorr, Fcorr] = myTimePowerSpectrum(ROItrace_corr, Fs);
h=plot(Fb, 10*log10(Pxxb),'color',[1 .5 0]); hold on;
plot(Fp, 10*log10(Pxxp),'r');
plot(Fcorr, 10*log10(Pxxcorr),'b');
xlabel('freq (hz)');
ylabel('power');
axis tight;
grid minor
marginplot;

subplot(224);
imagesc(imresize(ScaleFactor, [size(Up,1), size(Up,2)]));
hold on;
% plot(xidx/pixSpace, yidx/pixSpace,'wo');
contour(roi(:,:,1),[.5 .5],'w');
caxis([-1 1]*max(abs(ScaleFactor(:))));
colormap(colormap_blueblackred);
axis equal tight
colorbar
title('Correction Scale factor');

screen2png(fullfile(savePath,['HemoCorrect_sanityCheckFigs_' num2str(round(hemoFreq)) ' Hz']));%19/11/20
close(fig);

% % Vcorrfilt = detrendAndFilt(Vcorr, Fs); %5/1/18
% % svdViewer(Ub, DSb.Sv, Vcorrfilt, Fs) %blue  after filtering and detrending

writeUVtoNPY([], Vcorr, [], fullfile(savePath, ['svdTemporalComponents_corr']));
writeNPY(tb,  fullfile(savePath, 'svdTemporalComponents_corr.timestamps.npy'));
end


function [Pxx, F] = myTimePowerSpectrum(V, Fs)
L = length(V);
NFFT = 2^nextpow2(L);
[Pxx,F] = pwelch(V,[],[],NFFT,Fs);
end

