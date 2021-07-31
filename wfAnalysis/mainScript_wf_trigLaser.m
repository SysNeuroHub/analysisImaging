%% this is a script to do the following:
% detect onset of each stimulation in thorsync time using photodiode signal
% detect wf frame time in thorsync time
% obtain stimulus triggered time courses
% compute Direction selectivity index of each cell

%TODO
% analysis for ORsequence
% make a prompter to select expt info through GUI (use info in market
% server)

addpath(genpath('C:\npy-matlab'));
addpath(genpath('C:\Users\dshi0006\npy-matlab'));
addpath(genpath('C:\Users\Analysis\npy-matlab'));
addpath(genpath('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\master\Analysis\wfAnalysis\widefield'));
addpath(genpath('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\master\visbox\rigbox'));
addpath('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\master\Analysis');
addpath('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\master\Analysis\wfAnalysis');
addpath('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\Matteobox');

%% experiment
%ppbox notation
expt.subject = '97274';
expt.expDate = '2020-12-07_1';
expt.expNum = 1:;
bklightCtrl = 0;


%% SVD
nSV = 1000;
params.movieSuffix = 'blue';% 'purple''corr_dFF'; %'blue'  %U for corr_dFF can have NANs..
params.useCorrected = 1;

%% analysis
%highpassCutoff = 0.01; %[Hz]
%lowpassCutoff = []; %[Hz]
marginT = 3; %[s]
resizeS = 1;%0.25; %spatial rescaling factor
pixies = [150 92; 82 62; 122 24];%82 140; 105 145; 73 154; 98 145];
%[285 167; 295 139; 95 150; 55 137; 92 186; 113 162; 92 186; 188 144; 155 169;...
%     169 174; 188 144; 188 244; 194 235; 216 191; 263 150; ];


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

%movieWithTracesSVD(U,V,t)
%svdViewer(U, Sv, V, 35);
%
%% load mpep data
p = ProtocolLoad_wf(expt.subject,expt.expDate,expt.expNum); %3/6/20
figname = [dat.constructExpRef(expt.subject, thisDate, thisSeries, expt.expNum) '_'...
    p.xfile(1:end-2) '_' params.movieSuffix];

if contains(p.xfile, 'stimORsequence')
    doSequence = 1;
else
    doSequence = 0;
end

if ~doSequence
    stimSequence = getStimSequence(expt.subject, expt.expDate, expt.expNum);
    %stimSequence.seq
    %stimSequence.labels
    %stimSequence.paramValues
else
    [stimOnframes, stimOR] = grabORSequenceFrames( expt );
    orientations = unique([stimOR{1:3}]);
    
    %% under construction
    % %     stimSequence.labels=cell(p.nstim, 1);
    % %     stimSequence.seq=zeros(p.nrepeats*p.nstim, 1);
    % %
    % %     for ior = 1:length(orientations)
    % %         stimSequence.seq = find(stimOR{istim} == orientations(ior));
    % %     end
    % %     stimSequence.paramValues
end

if ~isfield(p,'pfiledurs')
    p.pfiledurs = p.pars(1,:)/10;
end

if ~doSequence
    respWin = [-marginT min(p.pfiledurs)+marginT]; %31/3/20
else
    tfparidx=find(strcmp(p.parnames, 'tf'));
    tf = unique(p.pars(tfparidx,:)/10); %[Hz]
    respWin = [-marginT 1/tf + marginT]; %12/11/20
end


%% detect onset of each stimulation in Timeline time using photodiode signal
if ~doSequence
    getAllPhdFlips = 0;
else
    getAllPhdFlips = [p.nstim p.nrepeats];
end
% expt = grabStimTimesWF(expt, getAllPhdFlips, [], [], bklightCtrl);
% % expt.stimTimes.onset
% % expt.stimTimes.offset
% % expt.stimTimes.onset
% % expt.stimTimes.frameTimes
AOTh = 0.1; %hack
expt = grabLaserTimesWF(expt,[],[],[],AOTh);

expt.stimTimes.onset = expt.laserTimes.onset;
expt.stimTimes.offset = expt.laserTimes.offset;
onDur = nanmedian(expt.laserTimes.offset - expt.laserTimes.onset);
%if there is no AO, substitute by trial onset
nolsrtrials = find(isnan(expt.stimTimes.onset));
expt.stimTimes.onset(nolsrtrials) = expt.laserTimes.trialOn(nolsrtrials);
expt.stimTimes.offset(nolsrtrials) = expt.stimTimes.onset(nolsrtrials) + onDur;

%check number of stimulation
if length(expt.stimTimes.onset) ~= p.nstim * p.nrepeats
    error(['Detected #stim onset:' num2str(expt.stimTimes.onset') '/ total #stim' num2str(p.nstim * p.nrepeats)]);
end

%% cf for 2p:
%[ expt ] = grabexpt.stimTimes.onset( expt, false, TSDir, thPhd, .9*ISI); %for 2p
% expt.expt.stimTimes.onset(expt.expNum).onset
% expt.expt.stimTimes.onset(expt.expNum).offset
% expt.stimTimes.onset = expt.expt.stimTimes.onset.onset; %[#stim onset x 1] time in thorsync time


if ~doSequence
    
    %% stimulus triggered movie
    pixelTuningCurveViewerSVD(U, V, t, expt.stimTimes.onset, stimSequence.seq, respWin);
    [avgPeriEventV, winSamps, periEventV] = ...
        eventLockedAvg(V, t, expt.stimTimes.onset, stimSequence.seq, respWin);
    %avgPeriEventV: icond x nSV x time
    %periEventV: event x nSV x time
    
    
    
    %% time-avg response & stimulus preference map
    preIdx = find(winSamps<0);
    postIdx = intersect(find(winSamps>0), find(winSamps < 0.5));%min(p.pfiledurs)));
    
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
    
    screen2png(fullfile(resultSaveDir,figname),gcf);
    close;
    
    
    %% single trial traces
    if ~isempty(pixies)
        for ff = 1:size(pixies,1)
            pixelSingleTrialTracesViewerSVD(U, V, t, mimg, p, expt.stimTimes, stimSequence, ...
                respWin, [pixies(ff,1) pixies(ff,2)],1);
            set(gcf, 'position',[0 0 1920 1080]);
            saveas(gcf, fullfile(resultSaveDir, [figname 'single_at' ...
                num2str([pixies(ff,1) pixies(ff,2)]) '.png']));
            close(gcf);
        end
    end

elseif doSequence
    tavgResp = nan(size(U,1), size(U,2), length(orientations), p.nstim);
    for istim = 1:p.nstim
        
        [~, stimCondIdx] = intersect(orientations, stimOR{istim});
        
        %         pixelTuningCurveViewerSVD(U, V, t, expt.stimTimes.frameTimes{istim}(stimOnframes{istim}), ...
        %             stimOR{istim}, respWin);
        
        % stimOnsetTimes = expt.stimTimes.frameTimes{istim}(stimOnframes{istim});
        stimOnsetTimes = expt.stimTimes.frameTimes{p.seqnums(istim)}(stimOnframes{istim}); %16/11/20
        stimIdentities = stimOR{istim};
        
        
        [avgPeriEventV, winSamps, periEventV] = ...
            eventLockedAvg(V,t, stimOnsetTimes, stimIdentities, respWin);
        
        %         [avgPeriEventV, winSamps, periEventV] = ...
        %             eventLockedAvg(V, t, expt.stimTimes.onset, stimSequence.seq, respWin);
        %         %avgPeriEventV: icond x nSV x time
        %         %periEventV: event x nSV x time
        
        
        %% time-avg response & stimulus preference map
        stimOnRatioIdx = find(strcmp(p.parnames, 'stimOnRatio'));
        stimOnRatio = unique(p.pars(stimOnRatioIdx,:)/100); %[0 1]
        stimOnDur = 1/tf * stimOnRatio;
        
        preIdx = find(winSamps<0);
        postIdx = intersect(find(winSamps>0), find(winSamps < stimOnDur)); %13/11/20
        
        %subtract by prestimulus in each condition
        tavgRespV = mean(avgPeriEventV(:,:,postIdx),3) - mean(avgPeriEventV(:,:,preIdx),3);
        %condition x nSV
        
        tavgResp(:,:,stimCondIdx,istim) = svdFrameReconstruct(U, tavgRespV');
        
        %% single trial traces
        
        %hack to run pixelSingleTrialTracesViwer
        p_c.nstim = length(orientations);%length(unique(stimIdentities));
        p_c.nrepeats = p.nrepeats;
        p_c.pfiledurs = stimOnDur*ones(p_c.nstim,1);
        stimTimes_c.onset = expt.stimTimes.frameTimes{p.seqnums(istim)}(stimOnframes{istim}); %16/11/20
        stimTimes_c.offset = stimTimes_c.onset + stimOnDur;
        for ii = 1:length(stimIdentities)
            stimSequence_c.seq(ii) = find(orientations == stimIdentities(ii)); %there must be an smarter way...
        end
        stimSequence_c.labels = num2cell(orientations);
        
        for ff = 1:size(pixies,1)
            pixelSingleTrialTracesViewerSVD(U, V, t, mimg, p_c, stimTimes_c, ...
                stimSequence_c, respWin, [pixies(ff,1) pixies(ff,2)], 0);
            set(gcf,'position',[0 0 1920 1080]);
            saveas(gcf, fullfile(resultSaveDir, [figname 'singleTrial istim'  ...
                num2str(istim) 'at ' ...
                num2str([pixies(ff,1) pixies(ff,2)]) '.png']));
            close(gcf);
        end
    end
    
    mtavgResp = squeeze(nanmean(nanmean(tavgResp,4),3));
    
    
    %% trial-avg across stim conditions
    panel = [];
    for istim = 1:p.nstim
        for istimcond = 1:length(orientations)
            panel(istimcond,istim) = subplot(p.nstim,length(orientations), istimcond + (istim-1)*length(orientations));
            imagesc(squeeze(tavgResp(:,:,istimcond,istim)) ...
                - mtavgResp);
            axis equal tight off;
            
            if istim == 1
                title(num2str(orientations(istimcond)));
            end
            if istimcond == 1
                ylabel(num2str(istim));
            end
        end
    end
    caxes(panel,[1 99],'indirect');mcolorbar;
    
    saveas(gcf, fullfile(resultSaveDir, [figname 'tavg.png']));
      
    
    %% orientation tuning
    oo = 2*pi/180 * orientations;
    tavgResp_eqprob = permute(squeeze(tavgResp(:,:,:,p.nstim)),[3 1 2]);
    [ pref, circvar, amp ] = circstats( tavgResp_eqprob, oo );
    % pref: [-pi pi]
    pref(pref<0) = pref(pref<0) + 2*pi;
    prefOR = 180/2/pi*pref;
    circvar(circvar<0)=0;
    
    figure('position',[0 0 900 900]);
    imagesc(prefOR,'alphadata',1-circvar);
    axis equal tight;
    grid minor;
    colormap(hsv);
    mcolorbar(gca,.5);
    saveas(gcf, fullfile(resultSaveDir, [figname 'prefOR.png']));
    close(gcf);
end
    
