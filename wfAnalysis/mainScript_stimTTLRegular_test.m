setPath_analysisImaging;

%% image used for DMD
 load('/home/daisuke/tmp/CCFBL_400x300pix_8x7grid/CCFBL_400x300pix_8x7grid_tmpC.mat',...
     'image4OI');


%% experiment
expt.subject = 'yamatotakeru';
expt.expDate = '2025-10-16_2';
expt.expNum = 2;
bklightCtrl = 0;

%% SVD
nSV = 1000;

%% analysis
marginT = .5; %[s]
resizeS = 1; %spatial rescaling factor

%% data location
thisDate = expt.expDate(1:10);
thisSeries = str2num(expt.expDate(12:end));
expPath = fileparts(dat.expFilePath(expt.subject, thisDate, thisSeries, expt.expNum, 'widefield','master'));
saveVpath = expPath;
mpepDir = dat.reposPath('main', 'master');

load(dat.expFilePath(expt.subject, thisDate, thisSeries, expt.expNum, 'Timeline', 'master'));


%% load wf data
params.movieSuffix = 'red';
params.useCorrected = 0;
tname = params.movieSuffix;


disp('Loading widefield data');
disp(expt)
[U, V, t, mimg, mask] = quickLoadUVt(expPath, nSV, saveVpath, params);
%[strobeOnTimes, strobeOffTimes, strobeDurs] = getStrobeTimes(Timeline, 'alloptrig');
 mask = imresize(mask, resizeS);
    mimg = imresize(mimg, resizeS);
    Ux = size(U,2);
    Uy = size(U,1);

    p = ProtocolLoad_wf(expt.subject,expt.expDate,expt.expNum); %3/6/20

    stimSequence = getStimSequence(expt.subject, expt.expDate, expt.expNum);
    %stimSequence.seq
    %stimSequence.labels
    %stimSequence.paramValues

    %expt = grabLaserTimesWF(expt); %NG for stimTTLRegular
expt = grabStimTTLRegularTimesWF(expt);
%DMD switch time
%laser amplitude in V
expt.DMDTimes.imageID = repmat(1:56,1,2); %HACK

    %% detect onset of each condition in Timeline time using photodiode signal
    %expt = grabStimTimesWF(expt, 0, [], [], bklightCtrl,0);
    

%% align by TTL within and across trials
respWin = [-0.2 0.6];
 [avgPeriEventV, winSamps, periEventV] = ...
        eventLockedAvg(V, t, expt.DMDTimes.onset, expt.DMDTimes.imageID, respWin);

%%
    preIdx = find(winSamps<0);
    postIdx = intersect(find(winSamps>0), find(winSamps < 0.4));%min(p.pfiledurs)));
    
    %subtract by prestimulus in each condition
    tavgRespV = mean(avgPeriEventV(:,:,postIdx),3) - mean(avgPeriEventV(:,:,preIdx),3);
    %condition x nSV
    
    tavgResp = svdFrameReconstruct(U, tavgRespV');
    %tavgResp = tavgResp - tavgResp(:,:,p.blankstims);%still blood vessel remains...
    
    % range_c = 100*squeeze(tavgResp./mimg);
    range_c = tavgResp;
    crange = prctile(range_c(:),[1 99]);
    
    nRows            = 7;%ceil(sqrt(p.nstim));
    nCols = 8;%ceil(p.nstim/nRows);
    figure('position',[0 0 1900 1200]);
    panel = [];
    ax = [];
    for iim = 1:max(expt.DMDTimes.imageID)
        panel(iim) = subplot(nRows,nCols,iim);
        imagesc(squeeze(tavgResp(:,:,iim)),'alphadata',mask);
        ax(iim) = gca;
        axis equal tight off;
        title(['imageID ' num2str(iim)]);
    end
    linkcaxes(ax, crange);

    linkaxes(panel);
    % screen2png(fullfile(resultSaveDir,'cameraImages'));
    screen2png('cameraImages');
    close;


%% compare recorded and stimulated images

tavgResp_c = tavgResp;
tavgResp_c(tavgResp > prctile(tavgResp(:),99)) = prctile(tavgResp(:),99);
iim = 10;
imshowpair(squeeze(tavgResp_c(:,:,iim)), squeeze(image4OI(:,:,iim)));
screen2png('camera_DMD_Images');
    