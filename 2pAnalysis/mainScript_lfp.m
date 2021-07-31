%% this is a script to do the following:
% detect onset of each stimulation in thorsync time using photodiode signal
% detect 2p frame time in thorsync time
% obtain stimulus triggered time courses
% compute Direction selectivity index of each cell

% make sure to add the paths
addpath('C:\Users\Analysis\Documents\MATLAB\ThorScanAnalysis');

%% experiment
%ppbox notation
expt.subject = 'c75_87';%'97061';
expt.expDate = '2020-04-05_3';%'2020-03-18';
expt.expNum = 2;


% expt.subject = 'c75_86';%'97061';
% expt.expDate = '2020-04-03_2';%'2020-03-18';
% expt.expNum  = 2;

%% analysis
marginT = 0.5; %[s]
fs_rs = 100; %down sampling rate[Hz]


%% data location
TSDir = 'D:\thorimagedata\'; % = ops0.RootStorage where ThorSync data is saved or copied
%S2PDir = 'D:\OutputSuite2p\';% = ops0.ResultsSavePath
mpepDir = 'Y:\'; %= ops0.mpepRootStorage


%% move LFP data to analysis PC
mpepInfo.animal = expt.subject;
mpepInfo.session = expt.expDate;
mpepInfo.experiment = expt.expNum;
[lookupTable_ti, lookupTable_ts] = renameThorImageFiles(mpepInfo, ...
    mpepDir, TSDir, 'Z:\ThorImage','Z:\ThorImage',1);


%% load mpep data
p = ProtocolLoad(expt.subject,expt.expDate,expt.expNum, 'donotload', mpepDir);
stimSequence = getStimSequence(expt.subject, expt.expDate, expt.expNum, mpepDir);
%stimSequence.seq
%stimSequence.labels
%stimSequence.paramValues

if ~isfield(p,'pfiledurs')
    p.pfiledurs = p.pars(1,:)/10;
end
respWin = [-marginT p.pfiledurs+marginT];

if mean(isnan(stimSequence.paramValues)) ~= 1
    % do nothing, default
else
    if strcmp(p.xfile, 'vmovie3sequentialGrating.x')% assuming xpose or ypose
            xValue = p.pars(strcmp(p.parnames, 'x'), :);
            yValue = p.pars(strcmp(p.parnames, 'y'), :);
        if numel(unique(xValue))>1 && numel(unique(yValue))==1
            stimSequence.paramValues = xValue;
        elseif numel(unique(xValue))==1 && numel(unique(yValue))>1
            stimSequence.paramValues = yValue;
        else 
            disp('needs coding!')
            disp('Press ENTER>>')
            pause
            stimSequence.paramValues=[];
        end
    else
            disp('needs coding!')
            disp('Press ENTER>>')
            pause
            stimSequence.paramValues=[];
    end
end

%% detect 2p frame time in thorsync time
[times_lfp, lfp] = grabLFP(expt, expt.expNum, 'all', TSDir);
sfreq = 1/median(diff(times_lfp));

%obtain powerspectrum of LFP
% [pxx,f] = pmtm(lfp, 3, length(lfp), sfreq);%slow
window = int16(2e4);
noverlap = int16(window/2);
frange = 0:500; %freqency range to show [hz]
[pxx,f] = pwelch(lfp-mean(lfp), window, noverlap, frange, sfreq);

figure('position',[0 0 500 900]);
subplot(311);
loglog(f,pxx);axis tight; grid minor
title([expt.subject ' ' expt.expDate ' ' num2str(expt.expNum)]);
xlabel('Frequency[Hz]');ylabel('lfp power');
% [~,fpassIdx] = max(pxx);
% fpass = [99 101];%[34 38];
% lfp_f = bandpass(lfp,fpass,sfreq); 
% lfp_p = (lfp_f - mean(lfp_f)).^2;


%% detect onset of each stimulation in thorsync time using photodiode signal
ISI = min(p.pfiledurs);
if isfield(p, 'minWait')
    ISI = ISI + min(p.minWait);
end
    
thPhd = .5;%'userInput';%
[ expt ] = grabStimTimes( expt, false, TSDir, thPhd, ISI*0.9);
% expt.stimTimes(expt.expNum).onset
% expt.stimTimes(expt.expNum).offset

%check number of stimulation
if  length(expt.stimTimes.onset) ~= p.nstim * p.nrepeats
    warning(['Detected #stim onset:' num2str(length(expt.stimTimes.onset')) '/ total #stim' num2str(p.nstim * p.nrepeats)]);
   
    stimOnCounts = length(expt.stimTimes.onset);
    stimOffCounts =  length(expt.stimTimes.offset);
    stimCounts = min(stimOnCounts, stimOffCounts);
    if stimOnCounts ~= stimOffCounts
        stimTimes = expt.stimTimes.onset(1:stimCounts);
    end
    
    % keep stimulus events when both photodiode and pfile are saved
    stimValidCounts = min(stimCounts, length(stimSequence.seq));
    stimTimes = stimTimes(1:stimValidCounts);
    stimSequence.seq = stimSequence.seq(1:stimValidCounts);
    if stimValidCounts < num2str(p.nstim * p.nrepeats)
        warning(['Only use ' num2str(stimValidCounts) 'events are used out of ' num2str(p.nstim * p.nrepeats)]);
    end
else
   stimTimes = expt.stimTimes.onset; %[#stim onset x 1] time in thorsync time 
end

stimTrace = zeros(length(times_lfp),1);
for ss = 1:length(stimTimes)
    tidx = find((times_lfp > expt.stimTimes.onset(ss)) .* (times_lfp < expt.stimTimes.offset(ss)));
    stimTrace(tidx) = 1;
end



%stackedplot(times_lfp, [lfp stimTrace]);
subplot(312);
for ss=1:size(stimTimes)
    rectangle('position',[expt.stimTimes.onset(ss) min(lfp) expt.stimTimes.offset(ss)-expt.stimTimes.onset(ss) max(lfp)-min(lfp)],...
        'edgecolor','none','facecolor',[.7 .7 .7]);
    hold on;
end
plot(times_lfp, lfp);
xlim([expt.stimTimes.onset(1)-0.5 expt.stimTimes.offset(p.nstim)+0.5]);%1st repeat
xlabel('Time from exp start [s]');
ylabel('lfp');

%% stimulus triggered time courses
[lfp_rs, times_rs] = resample(lfp, times_lfp, fs_rs);
[avgPeriEventV, winSamps, periEventV] = eventLockedAvg(lfp_rs', times_rs, stimTimes, stimSequence.seq, respWin);
avgPeriEventV = squeeze(avgPeriEventV);
% avgPeriEventV: nEventTypes x  nTimePoints
% winSamps: labels for the time axis, relative to the event times
eLabel = unique(stimSequence.seq); %1st dimension of avgPeriEventV
eLabel_paramValue = stimSequence.paramValues(eLabel);

% %subtract by prestimulus
preStimIdx = find(winSamps<0);
% avgPreStim = squeeze(nanmean(avgPeriEventV(:,preStimIdx),2));
% avgPeriEventV = avgPeriEventV - avgPreStim;
 
%visualization
subplot(313);
imagesc(winSamps, 1:length(eLabel), avgPeriEventV);
hold on;
line([0 0], [.5 length(eLabel)+.5],'linestyle','--','color',[.5 .5 .5]);
caxis(prctile(avgPeriEventV(:),[1 99]));
axis ij;
yticklabels(eLabel_paramValue);
xlabel('Time since stimulus onset[s]');
ylabel('stimulus');

figname = [expt.subject '_' expt.expDate '_' num2str(expt.expNum) '_summary.png'];
saveas(gcf,figname);
% close;

%% single trial traces
fig = figure('position',[ 680    68   1200   900]);
for istim = 1:p.nstim
    theseAxes(istim) = subplot(3,3,istim);
    
    theseEvents = find(stimSequence.seq == istim);
    
    avgPreStim = squeeze(nanmean(nanmean(periEventV(theseEvents,1,preStimIdx),1),3));
    plot(winSamps, squeeze(periEventV(theseEvents,1,:))-avgPreStim);hold on
    plot(winSamps, squeeze(avgPeriEventV(istim,:)-avgPreStim), 'k', 'linewidth',2);
     errorbar(winSamps, squeeze(avgPeriEventV(istim,:)-avgPreStim), ...
        squeeze(1/sqrt(length(theseEvents))*nanstd(periEventV(theseEvents,1,:))),...
        'k', 'linewidth',.5);
    title(eLabel_paramValue(istim));
    %ylim(prctile(periEventV(:),[.1 99.9]));
end
linkaxes(theseAxes);
xlabel('time since stimulus onset [s]');
saveas(fig, [expt.subject '_' expt.expDate '_' num2str(expt.expNum) '_singleTraces.png']);
% close(fig);

% fig = figure('position',[ 680    68   560   910]);
% for istim = 1:p.nstim
%     theseAxes(istim) = subplot(3,3,istim);
%     
%     theseEvents = find(stimSequence.seq == istim);
%     errorbar(winSamps, squeeze(avgPeriEventV(istim,:)), ...
%         squeeze(1/sqrt(length(theseEvents))*nanstd(periEventV(theseEvents,1,:))),...
%         'k', 'linewidth',.5);
%     title(eLabel_paramValue(istim));
%     ylim(prctile(periEventV(:),[.1 99.9]));
% end
% linkaxes(theseAxes);
% xlabel('time since stimulus onset [s]');
% saveas(fig, [expt.subject '_' expt.expDate '_' num2str(expt.expNum) '_avgTraces.png']);
% close(fig);

% stimIdx = eLabel(~isnan(eLabel_paramValue));
% blankIdx  = eLabel(isnan(eLabel_paramValue));
% stimValue = eLabel_paramValue(stimIdx);
% avgPeriEventV = avgPeriEventV(stimIdx, :, :);% - avgPeriEventV(blankIdx, :, :); %subtract by blank condition

