function [MyStack_ratioSum, ratioInfo, registerInfo, h, ttrace] = ...
    LoadOneStimulus_ratio(ServerDir, Cam, p, iStim, iTr, ResizeFactor, ratioInfo, registerInfo, varargin)
% make stackset of ratiometric. single stimulus, single repeat.
%
% S = StackSet.LoadOneStimulus_ratio(ServerDir, Cam, p, iStim, iTr, ResizeFactor)
% builds stackset of ratiometric, after noise-gain equalization for each pixel
% (Akemann 2012 JNP). It also returns a figure showing time trace of the ratiometric (using tools.showTtrace_Pspec) 
% Before running this function, Stackset for each FRET channel needs to be
% built
%
% Input: - ServerDir: Directory where StackSet of each camera is saved
%        - Cam: an Imager object. Cam.FileString will be overwritten
%        - p: p-file protocol;
%        - ResizeFactor: resize the image stacks according to ResizeFactor
%        - iStim: stimulus index. should be a scalar
%        - iTr: repeat number. should be a scalar
%
%  S = StackSet.LoadOneStimulus_ratio(ServerDir, Cam, p, iStim, iTr, ResizeFactor, ratioInfo)
%  specifies information for noise-gain equalization
%        - ratioInfo: information needed to process the ratiometric.
%           .firstFrameBase, .lastFrameBase: period to calculate noise-gain
%           correction factor. Empty > select from command line
%           .correctFreq: target frequency range [Hz] for noise-gain equalization.
%           Empty > [8 15]
%           .correctSpace: spatial filtering size for noise-gain
%           equalization. Empty > 15
%           .Fstop, .Fpass, Astop: paramters for high-pass filtering before
%           taking the ratio. If one of these are left empty, DO NOT apply
%           high-pass filtering
%
%  S = StackSet.LoadOneStimulus_ratio(ServerDir, Cam, p, iStim, iTr, ResizeFactor, ratioInfo, registerInfo)
%  applies image registration before taking the ratio.
%         - registerInfo: information needed to do image registration. See
%         makeRatioStack for details of ratioInfo/registerInfo
%           registerInfo.cam1, registerInfo.cam2: after image flipping
%
%  [~, ratioInfo] = StackSet.LoadOneStimulus_ratio(...)
%  returns ratioInfo, used to process noise-gain equalization. This option
%  is useful when ratioInfo is not inputted.
%
%  [~, ~, registerInfo] = StackSet.LoadOneStimulus_ratio(..., registerInfo)
%  returns registerInfo, used to process image registration.
%  registerInfo.Dx and Dy are added to the inputted registerInfo.
%  
%  [~, ~, ~, h] = StackSet.LoadOneStimulus_ratio(...)
%  returns a handle for the summary figure  
%
% See also: tools.Imager, fret.makeRatioStack, StackSet.LoadOneStimulus,  tools.showTtrace_Pspec
%
% 13-11-05 DS created
% 14-08-30 DS replace LoadStack with loadArr, now assume stackset of
% cam1&2 is already built and saved
% 15-01-27 DS added "suffix" as a varargin input

% TO DO:
% move cord to make a sanitiy-check figure to buildSingleTrStack
% move cord to temporaly align cam1 & 2 to makeratio
% modify suffix input to allow detrend and (high-pass) filtering
%
% 6/2/15 registerInfo.Dx, Dy are somehow not saved

%these should be recorded in somewhere
doSF = false;
ratioInfo.doHP = true;%default

% if ratioInfo.doHP
%     %% prepare high-pass filter
%     ratioInfo.Fstop = 0.5;  % Stopband Frequency
%     ratioInfo.Fpass = 1;  % Passband Frequency
%     ratioInfo.Astop = 5;%10;%14.1.4 for M131919 5;%13.12.26 for M131916 %20;   % Stopband Attenuation (dB)
%     %<these should be the input?
% end


    
if ~isempty(varargin)
    %registerInfo = varargin{1};
    
    vidx = [];
    for vv = 1 : length(varargin)
        if isstr(varargin{vv})
            vidx = [vidx vv];
        end
    end
    
    % any suffix for directory name. 27/1/15
    for vv = vidx
        if any(strfind(varargin{vv}, 'suffix'))
            suffix = varargin{vv+1};
            break
        end
    end

    if ~exist('suffix','var')
        suffix = '';
    end
end
    
if ~isequal(length(iStim), 1)
    error('iStim should be a scalar');
end

if isempty(ratioInfo)
    ratioInfo.correctFreq = [];
    ratioInfo.correctSpace = [];
    ratioInfo.firstFrameBase = [];
    ratioInfo.lastFrameBase = [];
    ratioInfo.rotateCam1 = true;
    ratioInfo.Fstop = [];
    ratioInfo.Fpass = [];
    ratioInfo.Astop = [];
end

% 13/5/14 DS
if isempty(ratioInfo.Fstop) || isempty(ratioInfo.Fpass) || isempty(ratioInfo.Astop)
    ratioInfo.doHP = false;
end

if ~isfield(ratioInfo,'correctFreq')
    ratioInfo.correctFreq = [];
end
if ~isfield(ratioInfo,'correctSpace')
    ratioInfo.correctSpace = [];
end
if ~isfield(ratioInfo,'firstFrameBase')
    ratioInfo.firstFrameBase = [];
end
if ~isfield(ratioInfo,'lastFrameBase')
    ratioInfo.lastFrameBase = [];
end

Cam.FileString = '_cam1';
% MyStack_cam1 = StackSet.LoadStacks( ServerDir, p, ResizeFactor, iStim, iTr, Cam);
% << should not use LoadStacks to avoid a loop, and remove ServerDir input
% MyStack_cam1 = StackSet.LoadOneStimulus(Cam, p, iStim, iTr, ResizeFactor, true);
StackDir = tools.getDirectory( ServerDir, p, ResizeFactor, iStim, iTr, Cam, suffix);
cam1name = fullfile(StackDir, sprintf('Stack%04d',iTr));

if exist([cam1name '.bin'], 'file') 
    [arr, MyStack_cam1] = tools.loadArr(cam1name); %29/8/2014 DS
    MyStack_cam1.Values = arr;
else
    error(['Could not find' cam1name '. Check if this stackset is already built']);
end
%re-applying fixmissingframes for interpolation
% mcam1 = squeeze(mean(mean(MyStack_cam1.Values)));
% nointerp_cam1 = find(~isnan(mcam1));
% ratioInfo.missStep_cam1 = find(isnan(mcam1));
% if ~isempty(ratioInfo.missStep_cam1)
%     MyStack_cam1.Values = ...
%         tools.fixMissingFramesFunc(MyStack_cam1.Values(:,:,nointerp_cam1), MyStack_cam1.TimeVec(nointerp_cam1), 'interp3');
% end

infoFileName = sprintf('infoStack%04d.mat', iTr);
if exist(fullfile(StackDir, infoFileName), 'file')
    load(fullfile(StackDir, infoFileName));
    interp_cam1 = info.timeInfo.interp;
    %[~, nointerp_cam1] = setxor(1:MyStack_cam1.nFrames, info.timeInfo.interp);
else
    display('Could not find info for missed frames');
    interp_cam1 = [];
    %nointerp_cam1 = 1:MyStack_cam1.nFrames;
end


Cam.FileString = '_cam2';
% MyStack_cam2 = StackSet.LoadStacks( ServerDir, p, ResizeFactor, iStim, iTr, Cam);
% MyStack_cam2 = StackSet.LoadOneStimulus(Cam, p, iStim, iTr, ResizeFactor, true);
StackDir = tools.getDirectory( ServerDir, p, ResizeFactor, iStim, iTr, Cam, suffix);
cam2name = fullfile(StackDir, sprintf('Stack%04d',iTr));
if exist([cam2name '.bin'], 'file') 
    [arr, MyStack_cam2] = tools.loadArr(fullfile(StackDir, sprintf('Stack%04d',iTr))); %29/8/14 DS
    MyStack_cam2.Values = arr;
else
   error(['Could not find' cam2name '. Check if this stackset is already built']);
end 
% mcam2 = squeeze(mean(mean(MyStack_cam2.Values)));
% nointerp_cam2 = find(~isnan(mcam2));
% ratioInfo.missStep_cam2 = find(isnan(mcam2));
% if ~isempty(ratioInfo.missStep_cam2)
%     MyStack_cam2.Values = ...
%         tools.fixMissingFramesFunc(MyStack_cam2.Values(:,:,nointerp_cam2), MyStack_cam2.TimeVec(nointerp_cam2), 'interp3');
% end

infoFileName = sprintf('infoStack%04d.mat', iTr);
if exist(fullfile(StackDir, infoFileName), 'file')
    load(fullfile(StackDir, infoFileName));
    if isfield(info, 'timeInfo')
        interp_cam2 = info.timeInfo.interp;
    else
        interp_cam2 = [];
    end
else
    display('Could not find info for missed frames');
    interp_cam2 = [];
end

if doSF
    ratioInfo.SFfilterSize = 3;
    disp(['(Master) Doing spatial filtering']);
    MyStack_cam1.Values = tools.averageFilter2(MyStack_cam1.Values, ratioInfo.SFfilterSize, 2);
    disp(['(Slave) Doing spatial filtering']);
    MyStack_cam2.Values = tools.averageFilter2(MyStack_cam2.Values, ratioInfo.SFfilterSize, 2);
end

%%

if ratioInfo.doHP
    Hd = fret.butterworthHP(MyStack_cam1.FrameRate, ratioInfo.Fstop, ratioInfo.Fpass, ratioInfo.Astop);
    
    disp(['(Master) Doing high-pass filtering at ' num2str(ratioInfo.Fstop) '[Hz]..']);
    tsteps1 = size(MyStack_cam1.Values,3);
    masterOffset = mean(MyStack_cam1.Values(:,:,ratioInfo.firstFrameBase:ratioInfo.lastFrameBase),3);
    [dum,MASTEROFFSET] = meshgrid(ones(tsteps1,1), masterOffset(:));
    MASTEROFFSET = reshape(MASTEROFFSET,size(masterOffset,1), size(masterOffset,2), tsteps1);
    
    MyStack_cam1.Values = filter(Hd, MyStack_cam1.Values - MASTEROFFSET, 3) + MASTEROFFSET;
    clear MASTEROFFSET
    
    disp(['(Slave) Doing high-pass filtering at ' num2str(ratioInfo.Fstop) '[Hz]..']);
    tsteps2 = size(MyStack_cam2.Values,3);
    slaveOffset = mean(MyStack_cam2.Values(:,:,ratioInfo.firstFrameBase:ratioInfo.lastFrameBase),3);
    [dum,SLAVEOFFSET] = meshgrid(ones(tsteps2,1), slaveOffset(:));
    SLAVEOFFSET = reshape(SLAVEOFFSET,size(slaveOffset,1), size(slaveOffset,2), tsteps2);
    
    MyStack_cam2.Values = filter(Hd, MyStack_cam2.Values - SLAVEOFFSET, 3) + SLAVEOFFSET;
    clear SLAVEOFFSET
end


%% rotation info %commented out on 14/4/14
% [ ratioInfo.t_concord, ratioInfo.fliplrCam1, ratioInfo.fliplrCam2, ratioInfo.flipudCam1, ratioInfo.flipudCam2] ...
%     = tools.LoadRotateInfo( ServerDir, p, ResizeFactor);


%for debugging
% mmaster2 = squeeze(mean(mean(MyStack_cam1.Values(110:130, 130:150,:))));
% mslave2 = squeeze(mean(mean(MyStack_cam2.Values(110:130, 130:150,:))));
% plot(mmaster2, 'r');hold on
% plot(mslave2 + mean(mmaster2-mslave2),'g');


%% align time stamps 
%get reference frame rate ...
% DataDir = fullfile(...
%     Cam.DataDir,[p.animal Cam.FileString],...
%     num2str(p.iseries), num2str(p.iexp) );
% [mFrameRate, mnFrames, mDur] = tools.getMedianTimes(DataDir, 'PCO');
% TimeStamp_align = (0:min(MyStack_cam1.nFrames, MyStack_cam2.nFrames)-1) / mFrameRate;


TimeStamp_align = MyStack_cam2.TimeVec; % 25/3/14 HACK ... this needs to be saved?

try
    [MyStack_cam1.Values, idx_cam1] = ...
        tools.alignTimeStamp(MyStack_cam1.Values, MyStack_cam1.TimeVec, TimeStamp_align, 'interp');
catch err %14/8/15 interp3 stucks when TimeVec is not strictly monotonic increasing
    [MyStack_cam1.Values, idx_cam1] = ...
        tools.alignTimeStamp(MyStack_cam1.Values, MyStack_cam1.TimeVec, TimeStamp_align, 'nearest');
end
MyStack_cam1.TimeVec = TimeStamp_align; %21/7/2015
MyStack_cam1.nFrames = length(idx_cam1); % 16/6/2014
[MyStack_cam2.Values, idx_cam2] = ...
    tools.alignTimeStamp(MyStack_cam2.Values, MyStack_cam2.TimeVec, TimeStamp_align);
MyStack_cam2.TimeVec = MyStack_cam2.TimeVec(idx_cam2); % 16/6/2014
MyStack_cam2.nFrames = length(idx_cam2); % 16/6/2014

[~, interp_cam1_aligned] = intersect(idx_cam1, interp_cam1);
[~, interp_cam2_aligned] = intersect(idx_cam2, interp_cam2);
if size(interp_cam1_aligned,1) > size(interp_cam1_aligned,2)
    interp_cam1_aligned = interp_cam1_aligned';
end
if size(interp_cam2_aligned,1) > size(interp_cam2_aligned,2)
    interp_cam2_aligned = interp_cam2_aligned';
end
ratioInfo.interp = unique([interp_cam1_aligned interp_cam2_aligned]); %interpolated frames (either cam1 or cam2) in steps


%% taking ratio
if exist('registerInfo','var') && ~isempty(registerInfo) %do image registration
    [MyStack_ratio, ratioInfo, registerInfo, MyStack_sum] = ...
        fret.makeRatioStack(MyStack_cam1, MyStack_cam2, ratioInfo, ...
        'registerInfo',registerInfo,'suffix', suffix);
    
    if ~isempty(registerInfo(2).registerImage)
        roix_trace = registerInfo(2).roix;
        roiy_trace = registerInfo(2).roiy;
    elseif ~isempty(registerInfo(1).registerImage)
        roix_trace = registerInfo(1).roix;
        roiy_trace = registerInfo(1).roiy;
    else
        error('cannot set roi for trace');
    end
else % without registration
    [MyStack_ratio, ratioInfo,~,MyStack_sum] = ...
        fret.makeRatioStack(MyStack_cam1, MyStack_cam2, ratioInfo,'suffix', suffix);
    
    roix_trace = round(MyStack_ratio.nCols*0.25) : round(MyStack_ratio.nCols*0.75);
    roiy_trace = round(MyStack_ratio.nRows*0.25) : round(MyStack_ratio.nRows*0.75);
end

MyStack_ratioSum = MyStack_ratio.AssignCondition(MyStack_sum.Values, 2);%2014/11/29

% commented out 19/5/14
% MyStack_ratio.TimeVec = TimeStamp_align;
% MyStack_ratio.nFrames = length(TimeStamp_align);

%% time traces for a sanity-check figure ... to be removed from here
roi_trace = zeros(MyStack_cam2.nRows, MyStack_cam2.nCols);
roi_trace(roiy_trace,roix_trace) = 1;

MyStack_cam1 = MyStack_cam1.Crop([roiy_trace(1) roiy_trace(end)],...
    [roix_trace(1) roix_trace(end)]);
if ~isempty(ratioInfo.masterBaseImage)
    masterBaseImage_crop = ratioInfo.masterBaseImage;
    masterBaseImage_crop([1:roiy_trace(1) roiy_trace(end):end],:) = [];
    masterBaseImage_crop(:,[1:roix_trace(1) roix_trace(end):end]) = [];
    MyStack_cam1.Values = tools.differentialSequence(MyStack_cam1.Values,...
        ratioInfo.firstFrameBase, ratioInfo.lastFrameBase, 1,MyStack_cam1.nFrames,...
        masterBaseImage_crop);
else
    MyStack_cam1.Values = tools.differentialSequence(MyStack_cam1.Values,...
        ratioInfo.firstFrameBase, ratioInfo.lastFrameBase);
end
MyStack_cam1 = MyStack_cam1.Trim([0.5 TimeStamp_align(end)-0.5]);
ttrace_cam1 = MyStack_cam1.SpaceAverages(1);

MyStack_cam2 = MyStack_cam2.Crop([roiy_trace(1) roiy_trace(end)],...
    [roix_trace(1) roix_trace(end)]);
if ~isempty(ratioInfo.slaveBaseImage)
    slaveBaseImage_crop = ratioInfo.slaveBaseImage;
    slaveBaseImage_crop([1:roiy_trace(1) roiy_trace(end):end],:) = [];
    slaveBaseImage_crop(:,[1:roix_trace(1) roix_trace(end):end]) = [];
    MyStack_cam2.Values = tools.differentialSequence(MyStack_cam2.Values,...
        ratioInfo.firstFrameBase, ratioInfo.lastFrameBase, 1,MyStack_cam2.nFrames,...
        slaveBaseImage_crop);
else
    MyStack_cam2.Values = tools.differentialSequence(MyStack_cam2.Values,...
        ratioInfo.firstFrameBase, ratioInfo.lastFrameBase);
end    
MyStack_cam2 = MyStack_cam2.Trim([0.5 TimeStamp_align(end)-0.5]);
ttrace_cam2 = MyStack_cam2.SpaceAverages(1);

MyStack_cam1.Values = [];
MyStack_cam2.Values = [];


%cf. LoadOneStimulusAverageRepeats
MyStack_ratioSum.Info.Stimulus = iStim;%otherwise Stimulus is named as "all"
MyStack_ratioSum.Info.RepeatList = iTr;
%MyStack_ratioSum.nConds = 1;%'Nrepeats'; % DS commented out on 27/8/2015

%MyStack_ratioSum.Description  = sprintf('%s-%d-%d Stimulus %d Average across repeats',...
%    p.animal,p.iseries,p.iexp, iStim);
MyStack_ratioSum.Description   = sprintf('%s-%d-%d-Stim %d RepeatList:%s',p.animal,...
    p.iseries,p.iexp,iStim, num2str(iTr));
%MyStack_ratioSum.SaveStacks(ServerDir, [], '_ratio');


%% traces of master, slave and ratio
MyStack_ratioCrop = MyStack_ratio.Crop([roiy_trace(1) roiy_trace(end)],...
    [roix_trace(1) roix_trace(end)]);
MyStack_ratioCrop = MyStack_ratioCrop.Trim([0.5 TimeStamp_align(end)-0.5]);
ttrace_time = MyStack_ratioCrop.TimeVec;
ttrace_ratio = MyStack_ratioCrop.SpaceAverages(1);

if any(isnan(ttrace_ratio)) || any(isinf(ttrace_ratio))
    warning('ratio trace includes nan or inf. Consider changing ROI.');
    ttrace_ratio(find(isnan(ttrace_ratio))) = 1;
    ttrace_ratio(find(isinf(ttrace_ratio))) = 1;
end

ttrace_length = min([length(ttrace_cam1) length(ttrace_cam2) length(ttrace_ratio)]);
ttrace(:,1) = ttrace_time(1:ttrace_length);
ttrace(:,2) = squeeze(ttrace_cam1(1:ttrace_length)) - nanmedian(squeeze(ttrace_cam1(1:ttrace_length)));%28/1/15
ttrace(:,3) = squeeze(ttrace_cam2(1:ttrace_length)) -  nanmedian(squeeze(ttrace_cam2(1:ttrace_length)));%28/1/15
ttrace(:,4) = squeeze(ttrace_ratio(1:ttrace_length)) -  nanmedian(squeeze(ttrace_ratio(1:ttrace_length)));%28/1/15

h1 = tools.showTtrace_Pspec(ttrace(:,2:4), ttrace(:,1), 'off');
subplot(3,2,1); 
title('Time courses');
ylabel('master dF/F');
subplot(3,2,3); 
ylabel('slave dF/F');
subplot(3,2,5); 
ylabel('ratio dR/R');
subplot(3,2,2);
title('Powerspectrum')

h2 = figure('visible','off');
subplot(2,2,1);
imagesc(ratioInfo.masterFactor);
axis equal tight;
caxis([1 2]);
pos = get(gca,'position');
colorbar;
hold on
contour(roi_trace,1,'color','w');
set(gca,'position',[pos(1) pos(2) pos(3) pos(4)]);
if isfield(ratioInfo, 'hbfreq')
    tname = sprintf('hb freq: %f [Hz]. Master factor',ratioInfo.hbFreq);
    title(tname);
end

if exist('registerInfo','var') && ~isempty(registerInfo)
    subplot(2,2,2);
    imagesc(registerInfo(2).registerImage);
    axis equal tight;
    pos = get(gca,'position');
    colorbar;
    hold on
    contour(roi_trace,1,'color','w');set(gca,'position',[pos(1) pos(2) pos(3) pos(4)]);
    title('register image cam2');
    
    subplot(2,2,4);
    [~,idxStart] = min(abs(MyStack_ratio.TimeVec - 0.5));
    [~,idxEnd] = min(abs(MyStack_ratio.TimeVec - (MyStack_ratio.TimeVec(end) -0.5)));
    plot( MyStack_ratio.TimeVec(idxStart:idxEnd), registerInfo(2).Dx(idxStart:idxEnd),...
        MyStack_ratio.TimeVec(idxStart:idxEnd), registerInfo(2).Dy(idxStart:idxEnd));
    yrange = prctile([registerInfo(2).Dx(idxStart:idxEnd) registerInfo(2).Dy(idxStart:idxEnd)], [1 99]);
    if ~yrange(1)==yrange(2) %14/12/1
        ylim(yrange);
    end
    xlim([MyStack_ratio.TimeVec(1) MyStack_ratio.TimeVec(end)]);
    legend('Dx','Dy');
    xlabel('Time [s]');
else
    registerInfo = [];
end

%from matteobox
h = mergefigs([h1; h2]);

scrsz = get(0,'ScreenSize');
set(h, 'position',[0 0 max(scrsz(3),1600) max(scrsz(4),1200)],'visible','off');

%close h1 h2

% Cam.FileString = '_ratio';
% StackDir = tools.getDirectory( ServerDir, p, ResizeFactor, iStim, iTr, Cam); 
% screen2png(fullfile(StackDir, ['TtracePspec_' num2str(iTr)]));
% save(fullfile(StackDir, ['TtracePspec_' num2str(iTr)]), 'ttrace');

end

