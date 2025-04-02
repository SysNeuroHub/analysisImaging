%% this is a script to do the following:
% load data from Market server
% detect onset of each stimulation in thorsync time using photodiode signal
% detect wf frame time in thorsync time
% obtain stimulus triggered time courses
% compute Direction selectivity index of each cell



setPath_analysisImaging;

%% experiment
expt.subject = 'izanami';
expt.expDate = '2025-01-30_1';
expt.expNum = 4;
bklightCtrl = 0;

%% SVD
nSV = 250;%1000;
params.movieSuffix = 'amber';% 'purple'
params.useCorrected = 1;

%% analysis
marginT = .5; %[s]
resizeS = 0.25; %spatial rescaling factor



%for estimation of preferred stim
n_boot = 10;%1 to see retinotopy only, 10 to compute VFS
use_method = 'max'; % max or com
screen_resize_scale = 3; %3 if max method


resultSaveDir = fullfile('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\AnalysisResult\wf',...
    [expt.subject '_' expt.expDate '_' num2str(expt.expNum)]);
mkdir(resultSaveDir);


%% data location
thisDate = expt.expDate(1:10);
thisSeries = str2num(expt.expDate(12:end));
expPath = fileparts(dat.expFilePath(expt.subject, thisDate, thisSeries, expt.expNum, 'widefield','master'));
saveVpath = expPath;
mpepDir = dat.reposPath('main', 'master');

load(dat.expFilePath(expt.subject, thisDate, thisSeries, expt.expNum, 'Timeline', 'master'));


%% load wf data
disp('Loading widefield data');
disp(expt)
[U, V, t, mimg, mask] = quickLoadUVt(expPath, nSV, saveVpath, params);
if isempty(t) %HACK when timeline strobes did not align with tiff frames
    [strobeOnTimes, strobeOffTimes, strobeDurs] = getStrobeTimes(Timeline, 'alloptrig');
    t = strobeOnTimes;
    t = t(1:size(V,2));
end
    
Fs = 1/median(diff(t));

U = imresize(U,resizeS);
mask = imresize(mask, resizeS);
mimg = imresize(mimg, resizeS);
Ux = size(U,2);
Uy = size(U,1);

%% dF/F NG
[newU, newV] = dffFromSVD(U, V, mimg);

%% temporal filtering V
fV = filtV(V, Fs, 0.01, 4);

%% load mpep data
p = ProtocolLoad_wf(expt.subject,expt.expDate,expt.expNum); %3/6/20
figname = [dat.constructExpRef(expt.subject, thisDate, thisSeries, expt.expNum) '_'...
    p.xfile(1:end-2) '_' params.movieSuffix];


if contains(p.xfile, 'stimORsequence')
    doSequence = 1;
else
    doSequence = 0;
end

stimSequence = getStimSequence(expt.subject, expt.expDate, expt.expNum);
%stimSequence.seq
%stimSequence.labels
%stimSequence.paramValues
if ~isfield(p,'pfiledurs')
    p.pfiledurs = p.pars(1,:)/10;
end

respWin = [-marginT min(p.pfiledurs)+marginT]; %31/3/20


%% detect onset of each stimulation in Timeline time using photodiode signal
expt = grabStimTimesWF(expt, 0, [], [], bklightCtrl,1);
% expt.stimTimes.onset
% expt.stimTimes.offset
% expt.stimTimes.frameTimes
%laserTh = 0.1;
%expt = grabLaserTimesWF(expt,[],[],[],laserTh);



%check number of stimulation
if length(expt.stimTimes.onset) ~= p.nstim * p.nrepeats
    error(['Detected #stim onset:' num2str(expt.stimTimes.onset') '/ total #stim' num2str(p.nstim * p.nrepeats)]);
end

%% hack for izanami 2025-01-30_2 exp4
lastValidTime = expt.stimTimes.offset(18*3);
[~, lastValidTidx] = min(abs(t - lastValidTime));

t = t(1:lastValidTidx);
V = V(:,1:lastValidTidx);
fV = fV(:,1:lastValidTidx);
newV = newV(:,1:lastValidTidx);

%% stimulus triggered movie
pixelTuningCurveViewerSVD(U, fV, t, expt.stimTimes.onset, stimSequence.seq, respWin);
[avgPeriEventV, winSamps, periEventV] = ...
    eventLockedAvg(fV, t, expt.stimTimes.onset, stimSequence.seq, respWin);
%avgPeriEventV: icond x nSV x time
%periEventV: event x nSV x time



%% time-avg response & stimulus preference map
preIdx = find(winSamps<0);
postIdx = intersect(find(winSamps>1), find(winSamps < 2));%min(p.pfiledurs)));

%subtract by prestimulus in each condition
tavgRespV = mean(avgPeriEventV(:,:,postIdx),3) - mean(avgPeriEventV(:,:,preIdx),3);
%condition x nSV

tavgResp = svdFrameReconstruct(U, tavgRespV');
%tavgResp = tavgResp - tavgResp(:,:,p.blankstims);%still blood vessel remains...

%% compute preferred stim in each pixel
response_grid = cell(1,p.nstim-1);
for istim = 1:p.nstim
    if istim == p.blankstims
        continue
    end
    %response_grid{istim} = tavgRespV(istim,:)';
    theseEvents = find(stimSequence.seq == istim);
    response_grid{istim} = reshape(mean(periEventV(theseEvents,:,postIdx),3) - mean(periEventV(theseEvents,:,preIdx),3), ...
        length(theseEvents),nSV)'; %[nSV %presentation]
end
prefMap = nanmean(prefStimSVD(U, response_grid,screen_resize_scale,n_boot,use_method),3);
prefMap(~mask)=nan;

nRows            = ceil(sqrt(p.nstim));
nCols = ceil(p.nstim/nRows);
figure('position',[0 0 1900 1200]);
panel = [];
for istim = 1:p.nstim-1
    panel(istim) = subplot(nRows,nCols,istim);
    imagesc(squeeze(tavgResp(:,:,istim)),'alphadata',mask);
    axis equal tight;
    title(stimSequence.labels{istim});
end
caxes(panel,[1 99],'indirect');mcolorbar;

panel_pref=subplot(nRows,nCols,p.nstim);

if contains(p.description, 'ori1')
    orientations = p.pars(p.activepars{1},1:p.nstim-1);
    oo = 2*pi/180 * orientations;
    tavgResp_permute = permute(squeeze(tavgResp(:,:,1:p.nstim-1)),[3 1 2]);
    [ pref, circvar, amp ] = circstats( tavgResp_permute, oo );
    %pref: [-pi pi]
    pref(pref<0) = pref(pref<0) + 2*pi;
    prefOR = 180/2/pi*pref;
    circvar(circvar<0)=0;
    
    imagesc(prefOR,'alphadata',1-circvar);
    axis equal tight;
    colormap(panel_pref,'hsv');
    mcolorbar(gca,.5);
    title('Preferred OR');
else
    imagesc(prefMap,'alphadata',mask);axis equal tight;
    title('Preferred stimulus');
    cb=colorbar;
    cb.Ticks = 1:p.nstim-1;
    cb.TickLabels = stimSequence.labels(cb.Ticks);
    colormap(panel_pref,'jet');
end

saveas(gcf,fullfile(resultSaveDir,figname),'fig');
close;


%% show movie
traces = prepareTimelineTraces(Timeline);
movieWithTracesSVD(U, fV, t, traces);

