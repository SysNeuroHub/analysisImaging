setPath_analysisImaging;

%% experiment

expt.subject = 'test';
expt.expDate = '2026-06-24_1';
expt.expNum = 1;
PCOFileName = 'dum.tif';

%  [frameTimes tl_flag] = getEyeFrameTimes(expt.subject, expt.expDate, expt.expNum);

%% eye camera images
info.expRef = dat.constructExpRef(expt.subject, expt.expDate(1:10), ...
    str2double(expt.expDate(12:end)), expt.expNum);
fullNames = dat.expFilePath(info.expRef, 'eyetracking');
[eyeFolder ,eyeFileStem, ~] = fileparts(fullNames{2});
 load(fullfile(eyeFolder, eyeFileStem));  %eyeLog

vReader = VideoReader(fullfile(eyeFolder, eyeLog.loggerInfo.Filename));
nFrames_eyecam = vReader.NumFrames;
%vidframes =  squeeze(read(vReader,[1 Inf])); %2220 frames

%% Timeline
thisDate = expt.expDate(1:10);
thisSeries = str2num(expt.expDate(12:end));
expPath = fileparts(dat.expFilePath(expt.subject, thisDate, thisSeries, expt.expNum, 'widefield','master'));
saveVpath = expPath;
mpepDir = dat.reposPath('main', 'master');

load(dat.expFilePath(expt.subject, thisDate, thisSeries, expt.expNum, 'Timeline', 'master'));

% ax=axes;
% for ii=1:15
%     ax(ii)=subplot(15,1,ii);
%     plot(tltime, Timeline.rawDAQData(:,ii));
%     grid on;
%     ylabel(Timeline.hw.inputs(ii).name)
%     if ii<15
%         axis off;
%     end
% end
% linkaxes(ax,'x');

tltime = Timeline.rawDAQTimestamps';
camStrobe_idx = strcmp({Timeline.hw.inputs.name}, 'camStrobe');
%         camStrobe_flicker = max(Timeline.rawDAQData(:,camStrobe_idx)) - ...
%             min(Timeline.rawDAQData(:,camStrobe_idx)) > 2;
camStrobe_th = max(Timeline.rawDAQData(:,camStrobe_idx))/2;
camStrobe_on = Timeline.rawDAQData(:,camStrobe_idx) < camStrobe_th; %[0 1]
camStrobeEv = trace2Event(camStrobe_on, tltime);
camStrobeTimes = mean(camStrobeEv,2);
nFrames_tl = numel(camStrobeTimes); %2219

%% load brain data
%  [nRows, nCols, TimeStamps, Stack] = LoadCustomPCO(FileName, 1, 1);
imstack = loadTiffStack(fullfile(eyeFolder,PCOFileName), 'tiffobj', 1);%2219
           
fprintf('%d frames in the eye cam video file\n', nFrames_eyecam);
fprintf('%d timestamps in the log file\n', length(eyeLog.TriggerData));
fprintf('%d frames in the TTL output from eye cam\n', nFrames_tl);
fprintf('%d frames in the brain cam video file\n', size(imstack,3));
 