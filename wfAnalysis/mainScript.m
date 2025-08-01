%% this is a script to do the following:
% load data from Market server
% detect onset of each stimulation in thorsync time using photodiode signal
% detect wf frame time in thorsync time
% obtain stimulus triggered time courses
% compute Direction selectivity index of each cell



setPath_analysisImaging;

%% experiment
<<<<<<< HEAD
expt.subject = 'achilles';
expt.expDate = '2025-05-21_1';
expt.expNum = 5;
bklightCtrl = 0;

%% SVD
nSV = 1000;%1000;
params.movieSuffix = 'amber';%'amber';% 'purple'
params.useCorrected = 0;
=======
expt.subject = 'test';
expt.expDate = '2025-06-18_2';
expt.expNum = 3;
bklightCtrl = 0;

%% SVD
nSV = 2000;%1000;
>>>>>>> origin/main

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
for icolor = 1:2
    
    switch icolor
        case 1
            params.movieSuffix = 'red';
            params.useCorrected = 0;
        case 2
            params.movieSuffix = 'amber';
            params.useCorrected = 0;
        case 3
            params.movieSuffix = 'amber';
            params.useCorrected = 1;
    end
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
    [newU, newV] = dffFromSVD(U, V, mimg);
    
    %% temporal filtering V
    fV = filtV(V, Fs, 0.01, 5);
    
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
%     pixelTuningCurveViewerSVD(U, V, t, expt.stimTimes.onset, stimSequence.seq, respWin);
    title(tname);
    
    

    %% time-avg response & stimulus preference map
    [avgPeriEventV, winSamps, periEventV] = ...
        eventLockedAvg(V, t, expt.stimTimes.onset, stimSequence.seq, respWin);
    
    preIdx = find(winSamps<0);
    [~,postIdx] = min(abs(winSamps - 0.65));%to examine laser artifact
    %      postIdx = intersect(find(winSamps>1), find(winSamps < 2));%min(p.pfiledurs)));
    
    %subtract by prestimulus in each condition
    tavgRespV = mean(avgPeriEventV(:,:,postIdx),3) - mean(avgPeriEventV(:,:,preIdx),3);
    %condition x nSV
    
    tavgResp = svdFrameReconstruct(U, tavgRespV');
    %tavgResp = tavgResp - tavgResp(:,:,p.blankstims);%still blood vessel remains...
    
    range_c = 100*squeeze(tavgResp./mimg);
    crange = max(abs(prctile(range_c(:),[1 99.99])));
    
    nRows            = 1;%ceil(sqrt(p.nstim));
    nCols = p.nstim;%ceil(p.nstim/nRows);
    figure('position',[0 0 1900 400]);
    panel = [];
    for istim = 1:p.nstim
        panel(istim) = subplot(nRows,nCols,istim);
        imagesc(100*squeeze(tavgResp(:,:,istim)./mimg),'alphadata',mask);
        axis equal tight;
        title(stimSequence.labels{istim});
        caxis([-crange crange]);
        [h,g]=mcolorbar;
    end
    linkaxes(panel);
%     xlim([40 110]); ylim([140 220])
    g.YLabel.String='dF/F [%]';
    screen2png(fullfile(resultSaveDir,figname));
    close;

end


% traces = prepareTimelineTraces(Timeline);
% movieWithTracesSVD(U, fV, t, traces);

% figure;
% imagesc(mimg);
% axis equal tight; colorbar;
% colormap(gray);
% caxis(prctile(mimg(:),[1 95]));
% title([expt.subject ' ' expt.expDate ' ' params.movieSuffix]);
