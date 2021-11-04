addpath('\\storage.erc.monash.edu\shares\R-MNHS-Syncitium\Shared\Daisuke\sandbox');
addpath(genpath('\\storage.erc.monash.edu\shares\R-MNHS-Syncitium\Shared\Daisuke\temp'));%stackset
if exist('C:\Users\dshi0006\git','dir')
    addpath(genpath('C:\Users\dshi0006\git'));
else
    addpath(genpath('C:\git'));    %dsbox,oi-tools,analysis-tools,analysisImaging, npy-matlab,neurostim, marmolab-stimuli
end

rootDir = '\\storage.erc.monash.edu\shares\R-MNHS-Syncitium\Shared\MarmosetData';
rootDir_OI = rootDir;%'E:\OI\';

%fixed parameters
batchSize = 1e3; %imaging loading


%% parameters for each exp
subject = 'rat1gou';
expDate = '20211026';
expName = '9'; %currently used only for imaging

ephysDirBase = fullfile(rootDir,subject,'oephys','experiment1/recording24');
stimName = [subject '.kalatsky.193725.mat'];

%% parameter for image analysis
rescaleFac = 0.25; %resize factor in x and y
polarity = -1; %whether to detect upward or downward deflection. -1 for IOS
%ephys digital only ... analog channel is NOT analyzed
trCh = 1;
camStrobeCh = 7;
%imaging preprocessing
cutoffFreq = 0.001; %[Hz] 0.02 is too high?
lpFreq = 2; %{hz] to reduce heart beat
%kalatsky analysis
filterFourierMaps = 0;
filtRF = 3;
delayPhase = 0;%-pi/2;%-3*pi/4;%pi; %added delay in radian
%0 for intrinsic imaging
%pi for calcium imaging
n_boot = 10;
xposRange = [30 120];%deg. to be obtained from cic record


imagingDirBase = fullfile(rootDir_OI,subject,'imaging');
saveDirBase = fullfile(rootDir,subject,'processed');
stimDirBase = fullfile(rootDir, subject, 'stimuli', expDate(1:4), ...
    expDate(5:6), expDate(7:8));

%no need to modify these
jsonFile = fullfile(ephysDirBase,'structure.oebin');
stimFile = fullfile(stimDirBase, stimName);


%% check if these data files exists
imagingDir_full = fullfile(imagingDirBase, expDate, expName);
saveDir_full = fullfile(saveDirBase, expDate, expName);
imageSaveName = fullfile(saveDir_full,...
    ['imageData_' expDate '_' expName '_resize' num2str(rescaleFac*100) '.mat']);

if ~exist(imagingDir_full, 'dir')
    error(['DONT EXIST imaging ' imagingDir_full]);
end
if ~exist(jsonFile, 'file')
    error(['DONT EXIST ephys ' jsonFile]);
end
if ~exist(stimFile, 'file')
    error(['DONT EXIST stim ' stimFile]);
end


%% load camera data ... assuming quickAnalysis is already done
disp(['Loading ' imageSaveName]);
load(imageSaveName,'imageData');
imageSize_r = size(imageData.imstack,[1 2]);

%% load stimuluda data
load(stimFile,'c');
contDir = get(c.gabor.prms.contDir, 'atTrialTime', inf);
duration = get(c.gabor.prms.duration,'atTrialTime',inf);
duration = duration(1)/1e3;
dur = get(c.gabor.prms.dur,'atTrialTime',inf,'trial',1);
blankDur = get(c.gabor.prms.bdur,'atTrialTime',inf,'trial',1);
stimDur = dur - blankDur;
calcWin = [0 duration];
stimLabels = contDir;
tgtFreq = 1/get(c.gabor.prms.dur); %[hz] for frequency analysis. Skip if empty
nrRepeats = c.blocks.nrRepeats;

%% load openEphys digital data
C = load_open_ephys_binary(jsonFile, 'continuous', 1);
%TODO: retrieve photodiode signal
fraxis = double(C.Timestamps); %this is frame number not time

%retrieve trial onset/offset
E = load_open_ephys_binary(jsonFile, 'events', 1);
srate = E.Header.sample_rate;
fr0 = fraxis(1);
taxis=1/srate*(fraxis-fr0); %[s]


theseEvents = find(E.ChannelIndex == trCh);
tr_ev =E.Data(theseEvents);
tr_fr = E.Timestamps(theseEvents)-fr0;

trON_fr = tr_fr(tr_ev>0); %frame number of trial onset
stimOnTimes = taxis(trON_fr); %time of trial offset
trOFF_fr = tr_fr(tr_ev<0); %frame number of trial onset
stimOffTimes = taxis(trOFF_fr); %time of trial offset


%nTrials = length(trON_fr);


%retrieve camera acquisition times
theseEvents = find(E.ChannelIndex == camStrobeCh);
cam_ev =E.Data(theseEvents);
cam_fr = E.Timestamps(theseEvents)-fr0;

camOnTimes = taxis(cam_fr);
Fcam = 1/median(diff(camOnTimes)); %camera effective sampling rate
disp(['#camera timestamps recorded in openEphys: ' num2str(numel(camOnTimes))]);

%% preprocessing
V = reshape(imageData.imstack, imageSize_r(1)*imageSize_r(2), []);%V: nPixels x nTimePoints
if ~isempty(cutoffFreq)
    meanV = mean(V,2);
    V =  hpFilt(V-meanV, Fcam, cutoffFreq); %cutoffFreq
    V = V + meanV;
end

if ~isempty(lpFreq)
    meanV = mean(V,2);%must be temporal avg
    V =  lpFilt(V-meanV, 1/median(diff(camOnTimes)), lpFreq); %cutoffFreq
    V = V + meanV;
end

%% trim data of each trial
[~, winSamps, periEventV, sortedLabels] ...
    = eventLockedAvg(V, camOnTimes, stimOnTimes, stimLabels, calcWin);
    
    
%% get AngleMaps of single trials
%< something wrong with ANGLEMAPS ... too much variability across trials
ANGLEMAPS = zeros(imageSize_r(1),imageSize_r(2),nrRepeats, 2);
ABSMAPS = ANGLEMAPS;
for iStim = 1:2
    
    theseEvents = find(stimLabels == stimLabels(iStim));
    
    for iTr = 1:nrRepeats
        thisEventV = squeeze(periEventV(theseEvents(iTr),:,:));%omit the 1st dimension
        
        Tensor = reshape(thisEventV,imageSize_r(1),imageSize_r(2),[]);

        %Tensor =  svdFrameReconstruct(U, thisEventV);
        %should subtract by mimg or else?
        %do detrend?
        
        Tensor(isnan(Tensor))=0;
        %Tensor = Tensor - mean(Tensor,3); %SHOULD NOT USE THIS - causing
        %more variability across trials
        %Tensor = Tensor - mean(mean(Tensor,1),2); %mean across pixels
        
        [ComplexMaps, AbsMaps, AngleMaps] = FourierMaps(Tensor, Fcam, tgtFreq);
        
        if filterFourierMaps
            realMaps_filt = wiener2(real(ComplexMaps), [filtRF filtRF]);
            imagMaps_filt = wiener2(imag(ComplexMaps), [filtRF filtRF]);
            
            ComplexMaps_filt = realMaps_filt + 1i * imagMaps_filt;
            AngleMaps = angle(ComplexMaps_filt);
        end
        %correction of angle produced by stackset.fouriermaps so that it
        %ranges between [-pi pi] (-pi = earliest, pi = latest)
        AngleMaps = -1*polarity*AngleMaps;
        
        %shift phase
        %use unwrap?
        %AngleMaps(AngleMaps < minPhase) = AngleMaps(AngleMaps < minPhase) + 2*pi;
        %AngleMaps = AngleMaps - 2*pi;
        
        AngleMaps = AngleMaps + delayPhase;
        AngleMaps(AngleMaps > pi) = AngleMaps(AngleMaps > pi) - 2*pi;
        AngleMaps(AngleMaps < -pi) = AngleMaps(AngleMaps < -pi) + 2*pi;
        
        ANGLEMAPS(:, :, iTr, iStim) = AngleMaps;
        ABSMAPS(:,:,iTr, iStim) = AbsMaps;
    end
end

%% combine opposite directions
randidx = randi(nrRepeats,2,n_boot);
xMap = (ANGLEMAPS(:,:,randidx(1,:),1) - ANGLEMAPS(:,:,randidx(2,:),2));
xMapm = nanmedian(xMap,3);
%xMapm(~mask)=nan;

mdelay_phase = squeeze(mean(mean(ANGLEMAPS+pi,4),3));
mdelay_sec = 1/tgtFreq/2/pi*mdelay_phase;

%% convert phase to visual angle
doublePhaseRange = [-2*pi 2*pi*stimDur/(blankDur+stimDur)];%in theory this should be used but
%doublePhaseRange = [-pi pi*stimDur/(blankDur+stimDur)];%in practice this must be correct
xMapm_vf = diff(xposRange)/diff(doublePhaseRange)*(xMapm - doublePhaseRange(1)) + xposRange(1); %[deg]

%% visualization
mmAbsMaps = squeeze(mean(mean(ABSMAPS,3),4));
for iStim = 1:2
    subplot(2,2,iStim);
    imagesc(nanmedian(ANGLEMAPS(:,:,:,iStim),3));%,'alphadata',mmAbsMaps/max(mmAbsMaps(:)))
    axis equal tight;
    caxis([-pi pi]);
    title(['phase cond' num2str(iStim) ' [rad]']);
end
mcolorbar;

subplot(223);
imagesc(mdelay_sec);
axis equal tight;
title('delay [s]')
mcolorbar;

subplot(2,2,4);
imagesc(xMapm_vf);%,'alphadata',mmAbsMaps/max(mmAbsMaps(:)));
axis equal tight;
title('visual angle [deg]')
mcolorbar;
