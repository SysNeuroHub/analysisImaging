%% this is a script to do the following:
% load data from Market server
% detect onset of each stimulation in thorsync time using photodiode signal
% detect wf frame time in thorsync time
% obtain stimulus triggered time courses
% compute Direction selectivity index of each cell



setPath_analysisImaging;

%% experiment
expt.subject = 'hercules';
expt.expDate = '2025-05-29_1';
expt.expNum = 3;
expt = grabScreenInfo(expt);
bklightCtrl = 0;

%% SVD
nSV = 1000;%1000;
params.movieSuffix = 'amber';
params.useCorrected =1;

%% analysis
marginT = .5; %[s]
resizeS = 0.25; %spatial rescaling factor
stimList = 1:3; %stimulus numbers to analyze
filterFourierMaps = 0; %filter fourier map in the complex domein
filtRF = 3;


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
fV = filtV(V, Fs, 0.01, 2);

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

%% stimulus triggered movie
pixelTuningCurveViewerSVD(U, V, t, expt.stimTimes.onset, stimSequence.seq, respWin);

%% frequency analysis (from getKalatskyRetinotopy.m)
pidx_tf = find(strcmp(p.parnames, 'tf1'));
tf = unique(p.pars(pidx_tf,stimList)) / 10; %[Hz]

dur = unique(p.pars(find(strcmp(p.parnames, 'dur'))))/10;%[s]

%get phase and amplitude at tf_corrected
[avgPeriEventV, winSamps, periEventV] = eventLockedAvg(V, t, expt.stimTimes.onset, ...
        stimSequence.seq, [0 dur]);
%avgPeriEventV: icond x nSV x time
%periEventV: event x nSV x time
    
ANGLEMAPS = zeros(Uy,Ux,p.nrepeats, p.nstim);
ABSMAPS = zeros(Uy,Ux,p.nrepeats, p.nstim);
for iStim = stimList
    
    theseEvents = find(stimSequence.seq == iStim);
    
    for iTr = 1:p.nrepeats
        thisEventV = squeeze(periEventV(theseEvents(iTr),:,:));%omit the 1st dimension
        
        
        Tensor =  svdFrameReconstruct(U, thisEventV);
        %should subtract by mimg or else?
        %do detrend?
        
        %Tensor = Tensor - mean(Tensor,3);
        Tensor = Tensor - svdFrameReconstruct(U,mean(mean(avgPeriEventV,1),3)');
        
        % temporal derivative
        Tensor = cat(3, diff(Tensor, 1, 3), zeros(Uy,Ux));

        [ComplexMaps, AbsMaps, AngleMaps] = FourierMaps(Tensor, Fs, tf(iStim));%2*tf(iStim));
        
        if filterFourierMaps
            realMaps_filt = wiener2(real(ComplexMaps), [filtRF filtRF]);
            imagMaps_filt = wiener2(imag(ComplexMaps), [filtRF filtRF]);
            
            ComplexMaps_filt = realMaps_filt + 1i * imagMaps_filt;
            AngleMaps = angle(ComplexMaps_filt);
        end
        %correction of angle produced by stackset.fouriermaps so that it
        %ranges between [-pi pi] (-pi = earliest, pi = latest)
        AngleMaps = -AngleMaps;
        %         AngleMaps(AngleMaps>2*pi) = AngleMaps(AngleMaps>2*pi) - 2*pi;
        %         AngleMaps(AngleMaps<0) = AngleMaps(AngleMaps<0) + 2*pi;
        
        %shift phase
        %use unwrap?
        %AngleMaps(AngleMaps < minPhase) = AngleMaps(AngleMaps < minPhase) + 2*pi;
        %AngleMaps = AngleMaps - 2*pi;
        
        AngleMaps(AngleMaps > pi) = AngleMaps(AngleMaps > pi) - 2*pi;
        AngleMaps(AngleMaps < -pi) = AngleMaps(AngleMaps < -pi) + 2*pi;
        
        ANGLEMAPS(:, :, iTr, iStim) = AngleMaps;
        ABSMAPS(:,:,iTr, iStim) = AbsMaps;
    end
end

figure;
for iStim = stimList
    ax(iStim) = subplot(2,3,iStim);
    imagesc(mean(ABSMAPS(:,:,:,iStim),3)); 
    caxis([0 0.25]);
    mcolorbar;
    title(['stim freq: ' num2str(tf(iStim))]);
    subplot(2,3,iStim+3);
    imagesc(angle(mean(exp(1i.*ANGLEMAPS(:,:,:,iStim)),3))); 
    %caxis([-pi/2 pi/2]); 
    mcolorbar;
end
screen2png(['freqAnalysis_dFdt_' expt.subject '_' expt.expDate '_' num2str(expt.expNum) '_' params.movieSuffix '_' num2str(params.useCorrected)]);
