%% this is a script to do the following:
% load data from Market server
% detect onset of each stimulation in thorsync time using photodiode signal
% detect wf frame time in thorsync time
% obtain stimulus triggered time courses
% compute Direction selectivity index of each cell



setPath_analysisImaging;

%% experiment

expt.subject = 'Confucious';
expt.expDate = '2026-05-21_1';
expt.expNum = 3;
bklightCtrl = 0;
polyScanPro = false;

params.movieSuffix = 'amber';
params.useCorrected = 0;

%% SVD
nSV = 2000;

%% analysis
marginT = 1; %[s]
resizeS = 0.1;%0.25; %spatial rescaling factor



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
tname = params.movieSuffix;
if params.useCorrected
    tname = [tname '-corrected'];
end

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
%[newU, newV] = dffFromSVD(U, V, mimg); %NG

%% temporal filtering V
fV = filtV(V, Fs, [], 6);

%% load mpep data
p = ProtocolLoad_wf(expt.subject,expt.expDate,expt.expNum); %3/6/20
figname = [dat.constructExpRef(expt.subject, thisDate, thisSeries, expt.expNum) '_'...
    p.xfile(1:end-2) '_' params.movieSuffix '_' num2str(params.useCorrected)];


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
expt = grabStimTimesWF(expt, 0, [], [], bklightCtrl,0);
laserTh = 1; %[V]
expt = grabLaserTimesWF(expt,[],laserTh);

%% DMD pattern index at the time of laser onset
nPtn = 2;%TODO: read from .seq file
[DMDIn_state, DMDOut_state] = getDMDState(Timeline, nPtn);


%% sanity check DMD pattern idx
[f_sanitycheck, NG_min] = DMD_ptn_sanityCheck(p.xfile, expt, Timeline, stimSequence, DMDIn_state, DMDOut_state);
screen2png([figname '_DMD_ptn_sanityCheck'], f_sanitycheck);


%check number of stimulation
if length(expt.stimTimes.onset) ~= p.nstim * p.nrepeats
    error(['Detected #stim onset:' num2str(expt.stimTimes.onset') '/ total #stim' num2str(p.nstim * p.nrepeats)]);
end

%% restore stimulation movie
% retrieve DMD pattern on OI space
load('M:\DMD images\whole_mode=4\whole_mode=4_Confucious.mat',...
    'image4OI');
load('M:\DMD images\whole_mode=4\whole_mode=4_Confucious.mat',...
    'image4OI_all_wCCF');

if ~polyScanPro
    image4OI = round(image4OI);
end
[DMD4OIsequence, tevent_DMD] = getDMD4OI_eventLockedAvg(image4OI, ...
    DMDIn_state, Timeline, t, stimSequence, expt.stimTimes, p.nstim, respWin);

%% stimulus triggered movie
% pixelTuningCurveViewerSVD(U, V, t, expt.laserTimes.onset, stimSequence.seq, respWin, 1);
% title(tname);


%% pixel-wise spectrogram
[avgPeriEventV, tevent, periEventV, sortLabels] = ...
    eventLockedAvg(fV, t, expt.stimTimes.onset, stimSequence.seq, respWin);

%% save stimulus triggered movie
for istim = 1:p.nstim
    clim = [-50 50];
    k2cmap = customcolormap([0 1],{'#00FFFF','#000000'});
    k2amap = customcolormap([0 1],{'#FFA500','#000000'});
    figure('position',[0 0 600 900]);
    ax(1) = subplot(211);
    ax(2) = subplot(212);
    inout = {squeeze(DMD4OIsequence(:,:,:,istim)),svdFrameReconstruct(U, squeeze(avgPeriEventV(istim,:,:)))};
    playMatrix(inout,{[0 2.5],[-50 50]},...
        ax,'colormap',{k2cmap,k2amap},...
        'framerate',60,...
        'alphadata',{1-image4OI_all_wCCF, imresize(1-image4OI_all_wCCF, [size(U,1) size(U,2)])},...
        'timevec',tevent_DMD,'saveVideo',true,...
        'mvName',[figname '_' num2str(istim)],...
        'quality',90,'bgBlack',false);
end

tic;
clear S_dB
for istim = 1:p.nstim
    [S_dB(:,:,:,:,istim),F,T] = spectrogramSVD(U,squeeze(avgPeriEventV(istim,:,:)),tevent, mask);
end
dur = toc

%% spectrogram, all pixels
iy = 1:size(U,1);
ix = 1:size(U,2);
figure('position',[0 100 800 500]);
for istim = 1:p.nstim
    ax(istim) = subplot(p.nstim,3,3*istim-2);
    imagescylog(T,F,squeeze(nanmean(nanmean(S_dB(iy,ix,:,:,istim),1),2)));
    title([num2str(p.pars(2,istim)) ' [mV]']);
end
hold on;
plot([0 p.pars(1,1)/10],[p.pars(3,1)/1e3 p.pars(4,1)/1e3],'k:'); text(p.pars(1,1)/10,p.pars(4,1)/1e3,'f');
plot([0 p.pars(1,1)/10],2*[p.pars(3,1)/1e3 p.pars(4,1)/1e3],'k:'); text(p.pars(1,1)/10,2*p.pars(4,1)/1e3,'2f');
plot([0 p.pars(1,1)/10],4*[p.pars(3,1)/1e3 p.pars(4,1)/1e3],'k:'); text(p.pars(1,1)/10,4*p.pars(4,1)/1e3,'4f');
ylim([F(2) F(end)]);
linkcaxes(ax, [-51 62]);
xlabel('time [s]'); ylabel('frequency [Hz]');
[~,h1] = mcolorbar;
h1.Label.String = 'dB';

%% normalised spectrogram, all pixels
mS_dB = repmat(squeeze(nanmean(nanmean(nanmean(nanmean(S_dB(iy,ix,:,:,:),1),2),4),5)),[size(T)]);

for istim = 1:p.nstim
    ax2(istim) = subplot(p.nstim,3,3*istim-1);
    imagescylog(T,F,squeeze(nanmean(nanmean(S_dB(iy,ix,:,:,istim),1),2))./mS_dB);
    caxis([0 2]);
end
xlabel('time [s]');
[~,h2]=mcolorbar;
h2.Label.String = 'dB/mean(dB)';

subplot(1,3,3);
imagesc(mimg);hold on; plot(ix,iy,'ro');colormap(gray);
hline(iy);vline(ix);
axis equal tight;


%% extract out the principal frequency of stimulation at each moment




