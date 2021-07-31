%% this is a script to do the following:
% detect onset of each stimulation in thorsync time using photodiode signal
% detect wf frame time in thorsync time
% obtain stimulus triggered time courses
% compute Direction selectivity index of each cell

%TODO
% analysis for ORsequence
% make a prompter to select expt info through GUI (use info in market
% server)

%addpath(genpath('C:\npy-matlab'));
%addpath(genpath('C:\Users\dshi0006\npy-matlab'));
%addpath(genpath('C:\Users\Analysis\npy-matlab'));
addpath(genpath('S:\biomed-physiology\imaging-data\MATLAB\functions\WideFieldAnalysis'));
addpath(genpath('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\master\Analysis\wfAnalysis\widefield'));
addpath(genpath('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\master\visbox\rigbox'));
addpath('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\master\Analysis');
addpath('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\master\Analysis\wfAnalysis');
addpath('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\Matteobox');

%% experiment
expt.subject = 'L4GCaMP6s_252';
expt.expDate = '2021-02-05_2';
expt.expNum = 4;


%% SVD
nSV = 1000;
params.movieSuffix = 'blue';% 'purple'
params.useCorrected = 1;

%% analysis
marginT = .5; %[s]
resizeS = .5;%0.25; %spatial rescaling factor



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


%% load wf data
disp('Loading widefield data');
disp(expt)
[U, V, t, mimg, mask] = quickLoadUVt(expPath, nSV, saveVpath, params);
Fs = 1/median(diff(t));

U = imresize(U,resizeS);
mask = imresize(mask, resizeS);
mimg = imresize(mimg, resizeS);
Ux = size(U,2);
Uy = size(U,1);

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
expt = grabStimTimesWF(expt, 0, [], [], 0,1);
% expt.stimTimes.onset
% expt.stimTimes.offset
% expt.stimTimes.onset
% expt.stimTimes.frameTimes
%laserTh = 0.1;
%expt = grabLaserTimesWF(expt,[],[],[],laserTh);

%check number of stimulation
if length(expt.stimTimes.onset) ~= p.nstim * p.nrepeats
    error(['Detected #stim onset:' num2str(expt.stimTimes.onset') '/ total #stim' num2str(p.nstim * p.nrepeats)]);
end


%% stimulus triggered movie
pixelTuningCurveViewerSVD(U, V, t, expt.stimTimes.onset, stimSequence.seq, respWin);
[avgPeriEventV, winSamps, periEventV] = ...
    eventLockedAvg(V, t, expt.stimTimes.onset, stimSequence.seq, respWin);
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




