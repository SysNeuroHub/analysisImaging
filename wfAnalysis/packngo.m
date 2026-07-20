function packed = packngo(expt, params)

% %% INPUT
% params.movieSuffix
% params.useCorrected
% params.nSV
% params.resizeS
%
% %% OUTPUT
% U, V, t, mimg, mask, expt
% DMD4OI
%
% expt.stimTimes
% expt.laserTimes
% expt.DMDstate
% expt.eyeTimes


%TODO: 
% add readme
% add function to restore movie from SVD


laserTh = 1; %[V]

%% data location
thisDate = expt.expDate(1:10);
thisSeries = str2num(expt.expDate(12:end));
expPath = fileparts(dat.expFilePath(expt.subject, thisDate, thisSeries, expt.expNum, 'widefield','master'));
saveVpath = expPath;
regDir = fullfile(fileparts(fileparts(expPath)), 'MR2DMDresult');
    
tldata = load(dat.expFilePath(expt.subject, thisDate, thisSeries, expt.expNum, 'Timeline', 'master'));
tltime = tldata.Timeline.rawDAQTimestamps';


%% CCF warped to OI 
imgPrefix = 'CCFBL_584x450pix';

allDirs = dat.paths;
imgDir = fullfile(fileparts(allDirs.mainRepository), 'DMD images', imgPrefix);

% load from the result of streo2DMD
if exist(fullfile(imgDir,[imgPrefix '_' expt.subject '.mat']), 'file')
    loaded = load(fullfile(imgDir,[imgPrefix '_' expt.subject '.mat']), 'image4OI_all_wCCF');
    ccf_oi = loaded.image4OI_all_wCCF;

else % OR create from scratch ... not yet working other than quetzal
    error('run stereo2DMD with imgPrefix=CCFBL_584x450pix');
    % load(fullfile(imgDir, [imgPrefix '_stereo']), 'camImg');
    % 
    % load(fullfile(regDir, 'Atlas_reg_info.mat'), 'mrangle','tform','tform2', 'image2','autoTform');
    % OIsize = size(image2);
    % 
    % ccf_stereo = imread(fullfile(imgDir, [imgPrefix '_stereo_all_wCCF.png']));
    % ccf_stereo(1,:) = 0; ccf_stereo(end,:) = 0; %HACK to remove image border
    % 
    % %requires antsenv and capability to run .sh
    % [~, ccf_oi] = applyStereo2DMD(double(ccf_stereo./max(ccf_stereo(:))), camImg.bregmapix, camImg.MmPerPixel, ...
    %     mrangle, tform, tform2, OIsize, camImg.MmPerPixel, regDir, autoTform);
end

packed.ccf = imresize(ccf_oi, params.resizeS);


%% load wf data
disp('Loading widefield data');
disp(expt)
[U, V, t, mimg, mask] = quickLoadUVt(expPath, params.nSV, saveVpath, params);
if isempty(t) %HACK when timeline strobes did not align with tiff frames
    [strobeOnTimes, strobeOffTimes, strobeDurs] = getStrobeTimes(tldata.Timeline, 'alloptrig');
    t = strobeOnTimes;
    switch param.movieSuffix
        case 'amber'
            t = t(1:2:end);
        case 'red'
            t = t(2:2:end);
    end
    t = t(1:size(V,2));
end

Fs = 1/median(diff(t));

[Uoriy, Uorix,~] = size(U);
packed.U = imresize(U, params.resizeS);
packed.V = V;
packed.t = t;
packed.mask = imresize(mask, params.resizeS);
packed.mimg = imresize(mimg, params.resizeS);

%% load mpep data
p = ProtocolLoad_wf(expt.subject,expt.expDate,expt.expNum); %3/6/20

expt.stimSequence = getStimSequence(expt.subject, expt.expDate, expt.expNum);
%stimSequence.seq
%stimSequence.labels
%stimSequence.paramValues
if ~isfield(p,'pfiledurs')
    p.pfiledurs = p.pars(1,:)/10;
end

if ismember(p.xfile, {'stimTTLChirp.x','stimTTLSwitch.x','stimTTLOsc.x'})
    doOptStim = 1;
else
    doOptStim = 0;
end


%% detect onset of each stimulation in Timeline time using photodiode signal
expt = grabStimTimesWF(expt, 0, [], [], params.bklightCtrl,0);


%% eye camera frame timinig
series = expt.expDate;
info.expRef = dat.constructExpRef(expt.subject, series(1:10), str2double(series(12:end)), expt.expNum);
fullNames = dat.expFilePath(info.expRef, 'eyetracking');
[eyeFolder ,eyeFileStem, ~] = fileparts(fullNames{2});
if exist(fullfile(eyeFolder, [eyeFileStem '.mat']) ,'file')>0
    expt.eyeTimes = grabEyeFrameTimesWF(expt.subject, expt.expDate, expt.expNum)';
end


if doOptStim
    expt = grabLaserTimesWF(expt,[],laserTh);
    %check number of stimulation
    if length(expt.stimTimes.onset) ~= p.nstim * p.nrepeats
        error(['Detected #stim onset:' num2str(expt.stimTimes.onset') '/ total #stim' num2str(p.nstim * p.nrepeats)]);
    end

    %% laser amplitude, loaded from Timeline
    laser_interp = getLaserAmp(tldata.Timeline);
    expt.laserAmp = interp1(tltime, laser_interp, t);

    %% DMD pattern index at the time of laser onset
    nPtn = p.pars(2);
    [DMDIn_state, DMDOut_state] = getDMDState(tldata.Timeline, nPtn);
    %% sanity check DMD pattern idx
    NG_min = DMD_ptn_sanityCheck(p.xfile, expt, tldata.Timeline, expt.stimSequence, DMDIn_state, DMDOut_state);
    expt.DMDidx = interp1(tltime, DMDOut_state, t,'nearest');
    if isempty(NG_min); disp('PASSED DMD_ptn_sanityCheck'); end

    %% retrieve DMD images uploaded by uploadSeqPng.m
    DMDfiles = dir(fullfile(expPath,'DMD', '*.png'));
    DMDfileNames = {DMDfiles.name};
    prefixNum = cellfun(@(x) sscanf(x, '%d', 1), DMDfileNames);
    % Sort by the numeric prefix
    [~, idx] = sort(prefixNum);
    DMDFileNames_sorted = DMDfileNames(idx);
    image4DMD = [];
    for ifile = 1:numel(DMDFileNames_sorted)
        image4DMD(:,:,ifile) =  imread(fullfile(expPath,'DMD',DMDFileNames_sorted{ifile}));
    end

    %% warp from DMD to OI space, using individual subject data stored in market
    load(fullfile(regDir, 'Atlas_reg_info.mat'), 'tform2');%, 'image2');

    %OIsize = size(image2);
    packed.DMDimg = imresize(DMD2OI(image4DMD, tform2, [Uoriy Uorix]), [size(packed.U,1), size(packed.U,2)]);

end


packed.expt = expt;



