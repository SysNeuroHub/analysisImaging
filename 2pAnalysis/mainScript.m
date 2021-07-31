%% this is a script to do the following:
% detect onset of each stimulation in thorsync time using photodiode signal
% detect 2p frame time in thorsync time
% obtain stimulus triggered time courses
% compute Direction selectivity index of each cell

% make sure to add the paths
addpath('C:\Users\Analysis\Documents\MATLAB\ThorScanAnalysis');

%% experiment
%ppbox notation
expt.subject = 'test';
expt.expDate = '2020-07-10_1';
expt.expNum = 2;
manualSave = 0; %if 1, manual saving mode. Only 1st repeat is analyzed

%% analysis
marginT = 1; %[s]
medFiltOrder = []; %if empty, dont apply median filter

resultSaveDir = [expt.subject '_' expt.expDate '_' num2str(expt.expNum)];
mkdir(resultSaveDir);


%% data location
TSDir = 'D:\thorimagedata\'; % = ops0.RootStorage where ThorSync data is saved or copied
S2PDir = 'D:\OutputSuite2p\';% = ops0.ResultsSavePath
mpepDir = 'Y:\'; %= ops0.mpepRootStorage
    

%% load mpep data
p = ProtocolLoad(expt.subject,expt.expDate,expt.expNum, 'donotload', mpepDir);
stimSequence = getStimSequence(expt.subject, expt.expDate, expt.expNum, mpepDir);
%stimSequence.seq
%stimSequence.labels
%stimSequence.paramValues
if ~isfield(p,'pfiledurs')
    p.pfiledurs = p.pars(1,:)/10;
end
if manualSave
    disp('Manual saving mode is selected. Only 1st repeat is analyzed.');
    p.nrepeats = 1;
    p.seqnums = p.seqnums(:,1);
    stimSequence.seq = stimSequence.seq(1:p.nstim);
end


respWin = [-marginT min(p.pfiledurs)+marginT]; %31/3/20

%% detect 2p frame time in thorsync time
[ frameTimes, TIdata ] = grabFrameTimes( expt, expt.expNum, 'all', TSDir, manualSave); %[1 x 2pframes]
nlayers = size(frameTimes,2);
%frameTimes [frame x layer]

frameTimes_rc = parseFrameTimes(frameTimes, p.nrepeats, p.blockEndDelay);
%frameTimes_rc{repeats}(frames, layer)

%% load 2p data
disp('Loading Suite2P data');
disp(expt)
%% load TI #frames ... probably theres quicker way ... load tiff files instead?
nFramesTI = [];
traces_r = cell(1,p.nrepeats);
for iplane = 1:nlayers 
    for irepeat = 1:p.nrepeats
        s2pname = sprintf('%s/%s/%s/%d/F_%s_%s_plane%d_repeat%d_proc.mat', S2PDir, ...
            expt.subject, expt.expDate, expt.expNum, expt.subject, expt.expDate, iplane, irepeat);
        load(s2pname);
        traces = dat.Fcell{iplane}';%[2pframes x cells] %F
        
        %% use only rois judged as cells
        cellIdx = find([dat.stat.iscell]); %only ROIs judged as cells
        %cellIdx = 1:length([dat.stat.iscell]); %all ROIs
        traces_r{irepeat} = traces(:,cellIdx);
        ncells = length(cellIdx);

         nFramesTI(iplane, irepeat) = size(traces_r{irepeat},1);
    end  
end

if manualSave
    nframes = min(size(frameTimes_rc{1},1), size(traces_r{1},1));
    frameTimes_r{1} = frameTimes_rc{1}(1:nframes,:);
    traces_r{1} = traces_r{1}(1:nframes,:);
else
    %% reduce thorsync frames to match thorimage frames
    frameTimes_r = adjustFrameTimes_repeat(frameTimes_rc, nFramesTI);
    %frameTimes_r{repeats}(frames, layer)
end


%% concatenate traces across repeats
iplane = 1; %not yet implemented
traces_a = [];
frameTimes_a = [];
for irepeat = 1:p.nrepeats
    traces_a = [traces_a; traces_r{iplane,irepeat}]; %[frames x cells]
    frameTimes_a = [frameTimes_a; frameTimes_r{iplane, irepeat}]; %[frames x 1]
end
%ncells(iplane) = size(traces_a,2);

%% median filter 30/3/20
if ~isempty(medFiltOrder)
    traces_a_cache = cat(1,flipud(traces_a),traces_a, flipud(traces_a));
    traces_a_cache = medfilt1(traces_a_cache, medFiltOrder);
    nframes_cache = size(traces_a,1);
    traces_a = traces_a_cache(nframes_cache+1:2*nframes_cache,:);
end

%% detect onset of each stimulation in thorsync time using photodiode signal
ISI = min(p.pfiledurs);
if isfield(p, 'minWait')
    ISI = ISI + min(p.minWait);
end
thPhd = 0.5;%'userInput';%
[ expt ] = grabStimTimes2p( expt, false, TSDir, thPhd, .9*ISI);
% expt.stimTimes(expt.expNum).onset
% expt.stimTimes(expt.expNum).offset

if manualSave
    expt.stimTimes.onset = expt.stimTimes.onset(1:p.nstim);
    expt.stimTimes.offset = expt.stimTimes.offset(1:p.nstim);
end
stimTimes = expt.stimTimes.onset; %[#stim onset x 1] time in thorsync time

%check number of stimulation
if length(stimTimes) ~= p.nstim * p.nrepeats
    error(['Detected #stim onset:' num2str(stimTimes') '/ total #stim' num2str(p.nstim * p.nrepeats)]);
end

%% traces of all cells
fig=figure;
tidx = find((frameTimes_a > expt.stimTimes.onset(1)-0.5) & ...
    (frameTimes_a < expt.stimTimes.offset(p.nstim)+0.5));
traces_a_stack = traces_a(tidx,:) - mean(traces_a(tidx,:));
distCells = max(max(abs(diff(traces_a_stack'))));-min(diff(traces_a_stack'));
traces_a_stack = traces_a_stack + (distCells*(1:ncells)); 
for ss=1:size(stimTimes)
    rectangle('position',[expt.stimTimes.onset(ss) min(traces_a_stack(:)) ...
        expt.stimTimes.offset(ss)-expt.stimTimes.onset(ss) max(traces_a_stack(:))-min(traces_a_stack(:))],...
        'edgecolor','none','facecolor',[.7 .7 .7]);
    hold on;
end
plot(frameTimes_a(tidx), traces_a_stack');
axis tight
yticks = median(traces_a_stack);
set(gca, 'ytick', yticks, 'yticklabel', 1:ncells);
xlabel('Time from exp start [s]');
ylabel('Cells');
saveas(fig, fullfile(resultSaveDir, 'timelapse.png')); 
close(fig);


%% stimulus triggered time courses
[avgPeriEventV, winSamps, periEventV] = eventLockedAvg(traces_a', frameTimes_a, stimTimes, stimSequence.seq, respWin);
% avgPeriEventV: nEventTypes x nCells x nTimePoints
% winSamps: labels for the time axis, relative to the event times
eLabel = unique(stimSequence.seq); %1st dimension of avgPeriEventV
eLabel_paramValue = stimSequence.paramValues(eLabel);

%27/3/20
% prestimIdx = find(winSamps<0);
% avgPeriEventV = avgPeriEventV - mean(avgPeriEventV(:,:,prestimIdx),3);

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
    fig = figure('position',[ 680    68   560   910]);
    for istim = 1:p.nstim
        theseAxes(istim) = subplot(3,3,istim);
        
        theseEvents = find(stimSequence.seq == istim);
        plot(winSamps, squeeze(periEventV(theseEvents,icell,:)));hold on
        plot(winSamps, squeeze(avgPeriEventV(istim,icell,:)), 'k', 'linewidth',2);
        title(eLabel_paramValue(istim));
    end
    linkaxes(theseAxes);
    xlabel('time since stimulus onset [s]');
    ylabel(num2str(icell));
    saveas(fig, fullfile(resultSaveDir, [num2str(icell) '.png']));
    close(fig);
end

%% FoV
roipos = [];
for r = 1:numel(dat.stat)
    roipos(r,1).ypos = dat.stat(r).ypix;% + min(dat.ops.yrange) - 1;
    roipos(r,1).xpos = dat.stat(r).xpix;% + min(dat.ops.xrange) - 1;
    roipos(r,1).iscell = dat.stat(r).iscell;
end
%% only show ROIs classified as cells
roipos = roipos(cellIdx);

I = dat.mimg_proc(:,:,2); %from redraw_meanimg

roidata = load_S2P_rois(roipos, I);
fig = figure;
roiOnMimg(roidata);
% mu = median(I(:));
% sd1 = mean(abs(I(I<mu+1e-7) - mu));
% sd2 = 1e-7 + mean(abs(I(I>mu-1e-7) - mu));
% caxis(mu + 5*[-sd1 sd2]);
saveas(gcf, fullfile(resultSaveDir, ['ROI.png']));
close(fig);

% cellmap = mapRoidata(roidata,0);
% imagesc(dat.ops.mimg);colormap(gray);hold on
% contour(cellmap,'color','r');
% axis equal tight

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
