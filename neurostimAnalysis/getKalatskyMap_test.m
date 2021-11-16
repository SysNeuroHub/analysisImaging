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
subject = 'rat2gou';%'rat1gou';
expDate = '20211104';%'20211026';
expName = '8';%'9'; %currently used only for imaging
ephysDirBase = fullfile(rootDir,subject,'oephys','experiment2/recording4');
stimName = [subject '.kalatsky.132143.mat'];

% ephysDirBase = fullfile(rootDir,subject,'oephys','experiment3/recording3');
% stimName = [subject '.kalatsky.180705.mat'];

%% parameter for image analysis
rescaleFac = 0.25; %resize factor in x and y
polarity = -1; %whether to detect upward or downward deflection. -1 for IOS

%for exp8
roi1=[42.3229166666667 84.6284722222222 10.6145833333333 11.0486111111111];
roi2=[63.3229166666667 84.4236111111111 11.0486111111111 11.4826388888889];
roi=cat(1,roi1,roi2);

%ephys digital only ... analog channel is NOT analyzed
trCh = 1;
camStrobeCh = 7;
%imaging preprocessing
cutoffFreq = 0.01; %[Hz] 0.02 is too high?
lpFreq = 2; %{hz] to reduce heart beat

%kalatsky analysis
filterFourierMaps = 1;
filtRF = 3;
%delayPhase = 0;%-pi/2;%-3*pi/4;%pi; %added delay in radian
%0 for intrinsic imaging
%pi for calcium imaging
n_boot = 10;
xposRange = [30 120];%deg. to be obtained from cic record

%TODO
%fix FourierMaps
%weigh angleMaps with trials of higher Amplitude

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


%% dF/F
V = (V - mean(V,2)) ./ mean(V,2);

%% trim data of each trial
[avgPeriEventV, winSamps, periEventV, sortedLabels, uniqueLabels] ...
    = eventLockedAvg(V, camOnTimes, stimOnTimes, stimLabels, calcWin);
AvgTensor = reshape(avgPeriEventV, size(avgPeriEventV,1), imageSize_r(1), ...
    imageSize_r(2),[]);
AvgTensor = AvgTensor - mean(AvgTensor,4);
mAngleMaps = [];
for iStim = 1:length(uniqueLabels)
    [~, ~, mAngleMaps(:,:,iStim)] = FourierMaps(squeeze(AvgTensor(iStim,:,:,:)), Fcam, tgtFreq);
end
delayPhase = mean(mAngleMaps,3);
delayPhase = 0.5*delayPhase;

%% get AngleMaps of single trials
%< something wrong with ANGLEMAPS ... too much variability across trials
ANGLEMAPS = zeros(imageSize_r(1),imageSize_r(2),nrRepeats, 2);
ABSMAPS = ANGLEMAPS;
for iStim = 1:length(uniqueLabels)
    
    theseEvents = find(stimLabels == uniqueLabels(iStim));
    
    for iTr = 1:nrRepeats
        thisEventV = squeeze(periEventV(theseEvents(iTr),:,:));%omit the 1st dimension
        
        %thisEventV = polarity * thisEventV; %for intrinsic imaging data
        %... NG due to bug? in FourirMaps
               
        Tensor = reshape(thisEventV,imageSize_r(1),imageSize_r(2),[]);
        
        mTensor = mean(mean(Tensor));
        nanFrames = find(isnan(mTensor));
        Tensor(:,:,nanFrames) = repmat(nanmean(Tensor,3),1,1,length(nanFrames));
        Tensor = Tensor - mean(Tensor,3); %SHOULD NOT USE THIS - causing
        %more variability across trials
        %Tensor = Tensor - mean(mean(Tensor,1),2); %mean across pixels
        
        %tmp
        S = StackSet;
        S.Values = Tensor;
        S.TimeVec = winSamps;
        S = S.FillDims;
        eventTimes = 0:1/tgtFreq:max(winSamps);
        tLims = [0 1/tgtFreq];
        Sm = S.MeanEvent(eventTimes, tLims);
        Sm.PlayCondition(1,prctile(Sm.Values(:),[1 99]),[],[],roi);
        
        [ComplexMaps, AbsMaps, AngleMaps] = FourierMaps(Tensor, Fcam, tgtFreq);
        %TODO: behvior of this function is weird:
        %         angle+ = fouriermaps(tesnsor)
        %         angle- = fouriermaps(-tensor)
        %         angle+ ~ -angle- + pi??
        
        if filterFourierMaps
            realMaps_filt = wiener2(real(ComplexMaps), [filtRF filtRF]);
            imagMaps_filt = wiener2(imag(ComplexMaps), [filtRF filtRF]);
            
            ComplexMaps_filt = realMaps_filt + 1i * imagMaps_filt;
            AngleMaps = angle(ComplexMaps_filt);
        end
        %correction of angle produced by stackset.fouriermaps so that it
        %ranges between [-pi pi] (-pi = earliest, pi = latest)
       
        
        AngleMaps = polarity * AngleMaps;
        
        %% shift phase
        AngleMaps = AngleMaps - delayPhase; %+?
        AngleMaps(AngleMaps > pi) = AngleMaps(AngleMaps > pi) - 2*pi;
        AngleMaps(AngleMaps < -pi) = AngleMaps(AngleMaps < -pi) + 2*pi;
        
        ANGLEMAPS(:, :, iTr, iStim) = AngleMaps;
        ABSMAPS(:,:,iTr, iStim) = AbsMaps;
        
    end
    
    images(ABSMAPS(:,:,:,iStim));
    screen2png(['singleAmp_cond' num2str(iStim)]);
    close(gcf);

    images(ANGLEMAPS(:,:,:,iStim));
    screen2png(['singleAngle_cond' num2str(iStim)]);
    close(gcf);
end

%% combine opposite directions
randidx = randi(nrRepeats,2,n_boot);
xMap = (ANGLEMAPS(:,:,randidx(1,:),2) - ANGLEMAPS(:,:,randidx(2,:),1));
xMapm = nanmedian(xMap,3);
%xMapm(~mask)=nan;

mdelay_phase = squeeze(mean(mean(ANGLEMAPS+pi+delayPhase,4),3));
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
