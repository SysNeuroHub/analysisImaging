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

bklightCtrl = 0;

tensors = zeros(175,176,90);
traces = [];
for aaa = 1:3
    switch aaa
        case 1
            expt.subject = 'L4GCaMP6s_250';
            expt.expDate = '2020-11-25_1';
            expt.expNum = 1;
        case 2
            expt.subject = 'L4GCaMP6s_250';
            expt.expDate = '2020-11-25_2';
            expt.expNum = 1;
        case 3
            expt.subject = 'L4GCaMP6s_250';
            expt.expDate = '2020-11-25_2';
            expt.expNum = 2;
    end
    
    %% SVD
    nSV = 1000;
    params.movieSuffix = 'blue';% 'purple''corr_dFF'; %'blue'  %U for corr_dFF can have NANs..
    params.useCorrected = 1;
    
    %% analysis
    %highpassCutoff = 0.01; %[Hz]
    %lowpassCutoff = []; %[Hz]
    marginT = 3; %[s]
    resizeS = 1;%0.25; %spatial rescaling factor
    %[285 167; 295 139; 95 150; 55 137; 92 186; 113 162; 92 186; 188 144; 155 169;...
    %     169 174; 188 144; 188 244; 194 235; 216 191; 263 150; ];
    
    
    %for estimation of preferred stim
    n_boot = 10;%1 to see retinotopy only, 10 to compute VFS
    use_method = 'max'; % max or com
    screen_resize_scale = 3; %3 if max method
    
    
    
    
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
    
    
    
    stimSequence = getStimSequence(expt.subject, expt.expDate, expt.expNum);
    %stimSequence.seq
    %stimSequence.labels
    %stimSequence.paramValues
    
    if ~isfield(p,'pfiledurs')
        p.pfiledurs = p.pars(1,:)/10;
    end
    
    respWin = [-0.5 min(p.pfiledurs)+marginT]; %31/3/20
    
    
    %% detect onset of each stimulation in Timeline time using photodiode signal
    getAllPhdFlips = 0;
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
    
    
    
    %% stimulus triggered movie
%     pixelTuningCurveViewerSVD(U, V, t, expt.stimTimes.onset, stimSequence.seq, ...
%         respWin,0);
    [avgPeriEventV, winSamps, periEventV] = ...
        eventLockedAvg(V, t, expt.stimTimes.onset, stimSequence.seq, respWin);
    %avgPeriEventV: icond x nSV x time
    %periEventV: event x nSV x time
    
    thiscond = 1;
    
    preIdx = find(winSamps<0);
    
    %thisV = squeeze(avgPeriEventV(thiscond,:,:)); %%NG!!!
    thisV = squeeze(mean(periEventV(stimSequence.seq==thiscond,:,:)));
    thisV = thisV - mean(thisV(:,preIdx),2);
    
    
    tensors(:,:,:,aaa) = svdFrameReconstruct(U,thisV);
    %     thisPixelU = squeeze(U(thisPixel(1),thisPixel(2),:));
    % traces(:,aaa)= thisPixelU'*squeeze(avgPeriEventV(thiscond,:,:));
    
    clear expt
end

%thisPixel = [88 88];%[106 66];%[91 73];%[100 78];%[73 62];%[66 66];


% traces = squeeze(tensors(pix(1),pix(2),:,:));
% ntraces = traces ./ max(traces,[],1);
%
% stackedplot(winSamps, ntraces)
tsize = size(tensors);
tensors2d_c = reshape(tensors,tsize(1)*tsize(2),tsize(3),tsize(4));
tensors2d = tensors2d_c(mask(:)==1,:,:);


for aaa = 1:3
    subplot(1,3,aaa);
    thisIm = squeeze(tensors2d(:,:,aaa));
    thisIm = thisIm / max(thisIm(:));
    imagesc(winSamps, 1:size(thisIm,1), thisIm);
    %caxis([-30 prctile(thisIm(:),[98])]);
    %caxis([prctile(thisIm(:),[5 98])]);
    caxis([-0.1 0.4]);
%     mcolorbar(gca,1,'southoutside');

end
colormap(1-gray);
screen2png('L4GCaMP6s_250_optogenResp_summary');

