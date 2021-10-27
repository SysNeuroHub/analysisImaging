function [MyStack_ratio, ratioInfo] = LoadOneStimulusAverageRepeats_ratio(ServerDir, Cam, p, ResizeFactor, iStim, ratioInfo)
% make ratio stackset. single stimulus, average repeats
%Input: - Cam, an Imager object. Cam.FileString will be overwritten
%        - p: p-file protocol;
%        - ResizeFactor: resize the image stacks according to ResizeFactor
%        - iStim: stimulus index. should be scalar
%        - ratioInfo: information needed to process the ratiometric.
%           .firstFrameBase, .lastFrameBase: period to calculate noise-gain
%           correction factor. Empty > select from command line
%           .correctFreq: target frequency range [Hz] for noise-gain equalization.
%           Empty > [8 15]
%           .correctSpace: spatialz` filtering size for noise-gain
%           equalization. Empty > 15
%           .
%Output:
%         - S: ratiometric signal after noise-gain equalization for each
%          pixel (Akemann et al. 2012 JNP)
%         - ratioInfo: information needed to process the ratiometric
%
% DS on 13.11.05
% function name is misleading because this mainly create and save the ratio

if ~isequal(length(iStim), 1)
    error('iStim should be a single number');
end

if isempty(ratioInfo)
    ratioInfo.correctFreq = [];
    ratioInfo.correctSpace = [];
    ratioInfo.firstFrameBase = [];
    ratioInfo.lastFrameBase = [];
    ratioInfo.rotateCam1 = true;
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

Cam.FileString = '_cam2';
MyStack_cam2 = StackSet.LoadStacks( ServerDir, p, ResizeFactor, iStim, [], Cam);

Cam.FileString = '_cam1';
MyStack_cam1 = StackSet.LoadStacks( ServerDir, p, ResizeFactor, iStim, [], Cam);


[ ratioInfo.t_concord, ratioInfo.fliplrCam1, ratioInfo.flipudCam1, ratioInfo.fliplrCam2, ratioInfo.flipudCam2] ...
    = tools.LoadRotateInfo( ServerDir, p, ResizeFactor, iStim, Cam);

%% prepare high-pass filter
elimstep = 50;
cutoffFreq = 0.5;
Fstop = 0.5;  % Stopband Frequency
Fpass = 1;  % Passband Frequency
Astop = 20;   % Stopband Attenuation (dB)
%<these should be the input?

Hd = fret.butterworthHP(MyStack_cam1.FrameRate, Fstop, Fpass, Astop);


disp(['(Master) Doing high-pass filtering at ' num2str(cutoffFreq) '[Hz]..']);
signal_cache = MyStack_cam1.Values;
nsx = size(signal_cache, 2);
nsy = size(signal_cache, 1);
for yy = 1:nsy;
    for xx = 1:nsx
        signal = squeeze(signal_cache(yy, xx, elimstep:end-elimstep));
        mergedsignal = double([flipud(signal); signal; flipud(signal)]) - mean(signal);
        mergedfiltered = filter(Hd,mergedsignal);
        MyStack_cam1.Values(yy, xx, :) = ...
            mergedfiltered(length(signal)+1-elimstep+1:2*length(signal)+elimstep) + mean(signal) - mean(mergedfiltered(length(signal)+1-elimstep+1:2*length(signal)+elimstep));
    end
end
clear signal_cache

disp(['(Slave) Doing high-pass filtering at ' num2str(cutoffFreq) '[Hz]..']);
signal_cache = MyStack_cam2.Values;
nsx = size(signal_cache, 2);
nsy = size(signal_cache, 1);
for yy = 1:nsy;
    for xx = 1:nsx
        signal = squeeze(signal_cache(yy,xx,elimstep:end-elimstep));
        mergedsignal = double([flipud(signal); signal; flipud(signal)]) - mean(signal);
        mergedfiltered = filter(Hd,mergedsignal);
        MyStack_cam2.Values(yy,xx,:) = ...
            mergedfiltered(length(signal)+1-elimstep+1:2*length(signal)+elimstep) + mean(signal) - mean(mergedfiltered(length(signal)+1-elimstep+1:2*length(signal)+elimstep));
    end
end
clear signal_cache

%for debugging
% mmaster = squeeze(mean(mean(MyStack_cam1.Values(150:200, 100:150,:))));
% mslave = squeeze(mean(mean(MyStack_cam2.Values(150:200, 100:150,:))));
% plot(mmaster, 'r');hold on
% plot(mslave + mean(mmaster-mslave),'g');

[MyStack_ratio, ratioInfo] = fret.makeRatioStack(MyStack_cam1, MyStack_cam2, ratioInfo);
%add ratioInfo.hbFreq, ratioInfo.masterFactor

MyStack_cam1.Values = [];
MyStack_cam2.Values = [];

%cf. LoadOneStimulusAverageRepeats
MyStack_ratio.Info.Stimulus = iStim;%otherwise Stimulus is named as "all"
MyStack_ratio.Info.RepeatList = 'Average';
MyStack_ratio.nConds = 1;
MyStack_ratio.Description  = sprintf('%s-%d-%d Stimulus %d Average across repeats',...
    p.animal,p.iseries,p.iexp, iStim);

%% save to the server
MyStack_ratio.SaveStacks(ServerDir, [], '_ratio');
%MyStack_ratio.Values = []; %oddly this line affects SaveStacks ..

%save info necessary for the ratiometric
AnimalDir   = fullfile(ServerDir, [MyStack_ratio.Info.animal '_ratio']);
SeriesDir   = fullfile(AnimalDir, sprintf('%03d', MyStack_ratio.Info.iseries));
ExpDir      = fullfile(SeriesDir, sprintf('%03d', MyStack_ratio.Info.iexp));
StackDir = fullfile(ExpDir, sprintf('Resize %d Stim %03d Repeat Average',...
    round(MyStack_ratio.ResizeFactor*100), MyStack_ratio.Info.Stimulus));
ThisFileName = sprintf('ratioInfoStack%04d.mat', MyStack_ratio.nConds);

save(fullfile(StackDir, ThisFileName), 'ratioInfo');
end

