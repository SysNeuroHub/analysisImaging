%% this is a script to do the following:
% detect onset of each stimulation in thorsync time using photodiode signal
% detect 2p frame time in thorsync time
% obtain stimulus triggered time courses
% compute Direction selectivity index of each cell

% previous version is saved in "mainScript_20200701" - HL, 20200701

% make sure to add the paths
addpath('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\master\Analysis');
addpath('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\master\Analysis\2pAnalysis');
addpath(genpath('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\master\visbox\rigbox'));

%% experiment
%ppbox notation
expt.subject = 'TIGRE2GCaMP6s_318';
expt.expDate = '2020-11-28_2';
expt.expNum = 3;
manualSave = 0; %if 1, manual saving mode. Only 1st repeat is analyzed

%% analysis
params.modality = 'Fcell';
params.useCells = 1;

marginT = 1; %[s]
medFiltOrder = 20; %if empty, dont apply median filter
subtractPre = 0; %15/12/20

savepng = 1; %if 1 save .png, else save .fig
resultSaveDir = fullfile('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\AnalysisResult\2p',...
    [expt.subject '_' expt.expDate '_' num2str(expt.expNum)]);
mkdir(resultSaveDir);


%% data location
TSDir = 'D:\thorimagedata\'; % = ops0.RootStorage where ThorSync data is saved or copied
S2PDir = 'D:\OutputSuite2p\';% = ops0.ResultsSavePath

%% copying data to the server 14/7/20
vaultDir = '\\vault-v2.erc.monash.edu.au\MNHS-dshi0006\Subjects'; %14/7/20
marketDir = 'M:\Subjects'; %14/7/20
TSDir_local = '\\ad.monash.edu\home\User006\dshi0006\Documents\TSDir_tmp';

tsnam_old = fullfile(vaultDir, expt.subject, expt.expDate, num2str(expt.expNum), 'Episode001.h5');
tsnam_new = fullfile(TSDir_local, expt.subject, expt.expDate, num2str(expt.expNum), 'Episode001.h5'); %12/1/21

if ~exist(tsnam_new,'file')
    mkdir(fullfile(TSDir_local, expt.subject, expt.expDate, num2str(expt.expNum)));
    copyfile(tsnam_old, tsnam_new);
    copyfile(fullfile(vaultDir, expt.subject, expt.expDate, num2str(expt.expNum), 'ThorRealTimeDataSettings.xml'),...
        fullfile(TSDir_local, expt.subject, expt.expDate, num2str(expt.expNum), 'ThorRealTimeDataSettings.xml'));
end

% migrateS2PtoServer(expt, TSDir, S2PDir, marketDir, vaultDir);

thisDate = expt.expDate(1:10);
thisSeries = str2num(expt.expDate(12:end));
figname = dat.constructExpRef(expt.subject, thisDate, thisSeries, expt.expNum);
if subtractPre %15/12/20
    figname = [figname '_subPre'];
end

%% load mpep data
p = ProtocolLoad_wf(expt.subject,expt.expDate,expt.expNum);%10/7/20 now this works if data is uploaded to Market!

if contains(p.xfile, 'stimORsequence') %12/1/21
    doSequence = 1;
else
    doSequence = 0;
end

if ~doSequence
    stimSequence = getStimSequence(expt.subject, expt.expDate, expt.expNum);%, mpepDir); %under Analysis not 2p
    %stimSequence.seq
    %stimSequence.labels
    %stimSequence.paramValues
else
    [stimOnframes, stimOR] = grabORSequenceFrames( expt );
    orientations = [0    15    30    45    60    75    90   105   120   135   150   165];

    %% under construction
    % %     stimSequence.labels=cell(p.nstim, 1);
    % %     stimSequence.seq=zeros(p.nrepeats*p.nstim, 1);
    % %
    % %     for ior = 1:length(orientations)
    % %         stimSequence.seq = find(stimOR{istim} == orientations(ior));
    % %     end
    % %     stimSequence.paramValues
end

    %% detect onset of each stimulation in thorsync time using photodiode signal
% for the data before 12/7/20
if datenum(expt.expDate(1:10)) <= datenum('2020-07-12')
    ISI = min(p.pfiledurs);
    if isfield(p, 'minWait')
        ISI = ISI + min(p.minWait);
    end
else
    % for the data after 12/7/20
    ISI = min(p.minWait);
end

if strcmp(expt.subject, '97274') && strcmp(expt.expDate, '2020-12-14_1') && (expt.expNum == 2)
    expt = grabStimTimes2ps( expt, true, TSDir_local, .8 ,5); %15/12/20    
else
    %thPhd = 0.5;
    %takes ~40sec
    if ~doSequence
        expt = grabStimTimes2p( expt, false, TSDir_local); %10/7/20
        % expt.stimTimes(expt.expNum).onset
        % expt.stimTimes(expt.expNum).offset
        
        %thr = 0.05; %20/11/20
        %expt = grabStimTimes2ps( expt, 0, marketDir, thr, ISI ); %20/11/20
    else
        expt = grabStimTimes2p( expt, true, TSDir_local); %21/8/20
    end
end

if manualSave
    expt.stimTimes.onset = expt.stimTimes.onset(1:p.nstim);
    expt.stimTimes.offset = expt.stimTimes.offset(1:p.nstim);
end
stimTimes = expt.stimTimes.onset; %[#stim onset x 1] time in thorsync time


%check number of stimulation
if length(stimTimes) ~= p.nstim * p.nrepeats
    error(['Detected #stim onset:' num2str(stimTimes') '/ total #stim' num2str(p.nstim * p.nrepeats)]);
end

if ~isfield(p,'pfiledurs')
    p.pfiledurs = p.pars(1,:)/10;
end
if manualSave
    disp('Manual saving mode is selected. Only 1st repeat is analyzed.');
    p.nrepeats = 1;
    p.seqnums = p.seqnums(:,1);
    stimSequence.seq = stimSequence.seq(1:p.nstim);
end


if ~doSequence
    respWin = [-marginT min(p.pfiledurs)+marginT]; %31/3/20
else
    respWin = [-marginT marginT]; %31/3/20
end

%takes ~40sec
[traces, expt.frameTimes, mimg, roi] = quickLoadS2Pt(expt, marketDir, p, params);

%% below for each plane
iplane = 1;
ncells = size(traces{iplane},1);

%% find F0 and dF/F0 16/12/20
tgtPrcTile = 20;
F0 = prctile(traces{iplane}, tgtPrcTile,2);
traces{iplane} = (traces{iplane}-F0)./F0;

%% median filter 30/3/20
if ~isempty(medFiltOrder)
    traces_cache = cat(2,fliplr(traces{iplane}),traces{iplane}, fliplr(traces{iplane}));
    traces_cache = medfilt1(traces_cache, medFiltOrder,[],2);
    nframes_cache = size(traces{iplane},2);
    traces{iplane} = traces_cache(:,nframes_cache+1:2*nframes_cache);
    clear traces_cache
end



if ~doSequence    
    %% traces of all cells of a specified repeat
    fig=figure('position',[0 0 800 1200]);
    irepeat = 10;
    theseStims = (1 + (irepeat-1)*p.nstim) : (irepeat*p.nstim);
    tidx = find((expt.frameTimes(iplane,:) > expt.stimTimes.onset(theseStims(1))-0.5) & ...
        (expt.frameTimes(iplane,:) < expt.stimTimes.offset(theseStims(end))+0.5));
    traces_stack = traces{iplane}(:,tidx) - mean(traces{iplane}(:,tidx));
    distCells = max(max(abs(diff(traces_stack))));%-min(diff(traces_stack'));
    traces_stack = traces_stack + (distCells*(1:ncells)');
    for ss = theseStims
        rectangle('position',[expt.stimTimes.onset(ss) min(traces_stack(:)) ...
            expt.stimTimes.offset(ss)-expt.stimTimes.onset(ss) max(traces_stack(:))-min(traces_stack(:))],...
            'edgecolor','none','facecolor',[.7 .7 .7]);
        hold on;
    end
    plot(expt.frameTimes(iplane,tidx), traces_stack');
    axis tight
    yticks = median(traces_stack');
    set(gca, 'ytick', yticks, 'yticklabel', 1:ncells);
    xlabel('Time from exp start [s]');
    ylabel('Cells');
    title(['repeat ' num2str(irepeat)]);
    if savepng
        saveas(fig, fullfile(resultSaveDir, 'timelapse.png'));
    else
        savefig(fig, fullfile(resultSaveDir, 'timelapse'));
    end
    close(fig);
    
    
    %% stimulus triggered time courses
    [avgPeriEventTrace, winSamps, periEventV] = eventLockedAvg(traces{iplane}, ...
        expt.frameTimes(iplane,:), stimTimes, stimSequence.seq, respWin);
    % avgPeriEventV: nEventTypes x nCells x nTimePoints
    % winSamps: labels for the time axis, relative to the event times
    eLabel = unique(stimSequence.seq); %1st dimension of avgPeriEventV
    eLabel_paramValue = stimSequence.paramValues(eLabel);
    
    if subtractPre
        prestimIdx = find(winSamps<0);
        avgPeriEventTrace = avgPeriEventTrace - mean(avgPeriEventTrace(:,:,prestimIdx),3);
        periEventV = periEventV - mean(periEventV(:,:,prestimIdx),3);
    end
    
    %visualization
    % for istim = 1:p.nstim
    %     subplot(3,3,istim)
    %     imagesc(winSamps, 1:ncells, squeeze(avgPeriEventV(istim,:,:)));
    %     caxis(prctile(avgPeriEventV(:),[1 99]));
    %     title(eLabel_paramValue(istim));
    % end
    % xlabel('time since stimulus onset');
    % ylabel('cell ID');
    % figname = [expt.subject '_' expt.expDate '_' num2str(expt.expNum) '.png'];
    % saveas(gcf,figname);
    
    
    for icell = 1:ncells
        fig = figure('position',[0 0 560 910]);
        for istim = 1:p.nstim
            theseAxes(istim) = subplot(3,3,istim);
            
            theseEvents = find(stimSequence.seq == istim);
            plot(winSamps, squeeze(periEventV(theseEvents,icell,:)));hold on
            plot(winSamps, squeeze(avgPeriEventTrace(istim,icell,:)), 'k', 'linewidth',2);
            title(eLabel_paramValue(istim));
        end
        linkaxes(theseAxes);
        xlabel('time since stimulus onset [s]');
        ylabel(num2str(icell));
        if savepng
            saveas(fig, fullfile(resultSaveDir, [figname '_' num2str(icell) '.png']));
        else
            savefig(fig, fullfile(resultSaveDir, [figname '_' num2str(icell)]));
        end
        close(fig);
        set(gcf,'Name',num2str(icell),'NumberTitle','off')
    end
    
elseif doSequence %currently only for ORsequence
  
   
%     U = reshape(eye(ncells),ncells,1,ncells); %dummy U to use prefStimSVD and svdFrameReconstruct
% 
%     [avgPeriEventV, winSamps, periEventV, sortedLabels] = ...
%         eventLockedAvgSVD(U, traces{iplane}, frame_t, ...
%         flip_times{eventIdx}(stimOnframes), stimOR, surround_window);
%     
%     periEventV = permute(periEventV, [2 3 1]);%[3 1 2]);
%     periEventMov = svdFrameReconstruct(Udf(:,:,1:nComps), periEventV(1:nComps,:));
%     periEventMov = reshape(periEventMov, [Ur Uc length(winSamps) length(stimOR)]);
%     
%     periEventTrace = squeeze(mean(mean(periEventMov(roi(1):roi(1)+roi(3),roi(2):roi(2)+roi(4),:,:))));
%     %time x trial

istim = p.nstim;%
    [avgPeriEventTrace, winSamps, periEventV] = ...
        eventLockedAvg(traces{iplane}, expt.frameTimes(iplane,:), ...
        expt.stimTimes.frameTimes{istim}(stimOnframes{istim}), stimOR{istim}, respWin);
    % avgPeriEventV: nEventTypes x nCells x nTimePoints

    %      eLabel = unique(stimSequence.seq); %1st dimension of avgPeriEventV
    %     eLabel_paramValue = stimSequence.paramValues(eLabel);

    nRows = ceil(sqrt(length(orientations)));
    nCols = ceil(length(orientations)/nRows);

    for icell = 1:ncells
        fig = figure('position',[0 0 560 910]);
        for istimSeq = 1:length(orientations)
            
            theseAxes(istimSeq) = subplot(nRows,nCols,istimSeq);
            
            %theseEvents = find(stimSequence.seq == istimSeq);
            theseEvents = find(stimOR{istim} == orientations(istimSeq));
            plot(winSamps, squeeze(periEventV(theseEvents,icell,:)));hold on
            plot(winSamps, squeeze(avgPeriEventTrace(istimSeq,icell,:)), 'k', 'linewidth',2);
            title(orientations(istimSeq));
        end
        linkaxes(theseAxes);
        xlabel('time since stimulus onset [s]');
        ylabel(num2str(icell));
        if savepng
            saveas(fig, fullfile(resultSaveDir, [num2str(icell) '.png']));
        else
            savefig(fig, fullfile(resultSaveDir, [num2str(icell)]));
        end
        close(fig);
        set(gcf,'Name',num2str(icell),'NumberTitle','off')
    end
    
end
%% FoV modified from roiOnMimg.m
fig=figure;
imagesc(mimg(:,:,iplane));colormap(gray);hold on;
contour(sum(roi{iplane},3),'color',[.9 .4 .4]);
axis equal tight

for icell = 1:ncells
    c = regionprops(roi{iplane}(:,:,icell),'centroid');
    text(c.Centroid(1),c.Centroid(2),num2str(icell),'color','g');
end
if savepng
    saveas(gcf, fullfile(resultSaveDir, [figname '_ROI.png']));
else
    savefig(gcf, fullfile(resultSaveDir, [figname '_ROI']));
end

close(fig);


%% average response across cells
fig = figure('position',[0 0 560 910]);
for istim = 1:p.nstim
    theseAxes(istim) = subplot(3,3,istim);
    
    theseEvents = find(stimSequence.seq == istim);
    plot(winSamps, squeeze(mean(periEventV(theseEvents,:,:),2)));hold on
    plot(winSamps, squeeze(median(mean(periEventV(theseEvents,:,:),2),1)), 'k', 'linewidth',2);
    %plot(winSamps, squeeze(mean(avgPeriEventTrace(istim,:,:),2)), 'k', 'linewidth',2);
    title(eLabel_paramValue(istim));
end
linkaxes(theseAxes);
xlabel('time since stimulus onset [s]');
ylabel(num2str(icell));
if savepng
    saveas(fig, fullfile(resultSaveDir, [figname '_avgCells.png']));
else
    savefig(fig, fullfile(resultSaveDir, [figname '_avgCells.png']));
end
close(fig);


% % figure;
% % ImageInfos = Marius2Takashi(dat);
% % nCell = size(ImageInfos.filterMatrix,1);
% % imagesc(dat.ops.mimg);
% % hold on;
% % for iCell = 1:nCell
% %     contour(squeeze(ImageInfos.filterMatrix(iCell,:,:)),'color','w');
% %     %         pause(0.1);
% %     %         contour(squeeze(thisImageInfos.filterMatrix_neuropil(iCell,:,:)),'color','g');
% %     %         plot(thisImageInfos.centerPosition(iCell,1), thisImageInfos.centerPosition(iCell,2), 'rx');
% % end
% % axis equal tight


% % %% DSI(n) of each cell
% % stimIdx = eLabel(~isnan(eLabel_paramValue));
% % blankIdx  = eLabel(isnan(eLabel_paramValue));
% % stimValue = eLabel_paramValue(stimIdx);
% % avgPeriEventV = avgPeriEventV(stimIdx, :, :);% - avgPeriEventV(blankIdx, :, :); %subtract by blank condition
% % %cf. http://tips.vhlab.org/data-analysis/measures-of-orientation-and-direction-selectivity
% % tavgPeriEventV = mean(avgPeriEventV, 3); %stimuli x cells
% %
% % [prefOriResp,prefOriIdx] = max(tavgPeriEventV, [],1); %cells x 1
% % prefOri = stimValue(prefOriIdx); % 1 x cells
% % orthOri = mod(prefOri+180, 360); % 1 x cells
% % orthOriIdx = zeros(1,ncells);  %cells x 1
% % orthOriResp = zeros(1,ncells,1);
% % for icell = 1:ncells
% %     [~, orthOriIdx(icell)] = min(abs(stimValue - orthOri(icell)));
% %     orthOriResp(icell) = tavgPeriEventV(orthOriIdx(icell),icell);
% % end
% %
% % DSIn = (prefOriResp - orthOriResp) ./ (prefOriResp + orthOriResp);
% %
% %
