function [ times, TIdata ] = grabFrameTimes( expt, expr, frames, TSDir, manualSave)
% [ times ] = grabFrameTimes( expt, expr, frames, TSDir)
% returns times of TTL signals going from LOW to HIGH
% for 2p frame time (default), return times of each plane as a cell (planes x 1)
% for other TSFieldNames ,return an array
%13/1/20 created from grabFrameTimes for Timeline version
% 10/2/20 deleted 5th input

subject = expt.subject;
expDate = expt.expDate;

TSFieldName = 'FrameOut';
if nargin < 5
    manualSave = 0;
end
if nargin < 4
    TSDir = 'D:\thorimagedata\'; %where ThorSync data is saved or copied
end

if nargin < 3
    frames = 'all';
end


%timelineFile = strcat(expDate,'_',num2str(expr),'_',subject,'_','Timeline.mat');
fnam = fullfile(TSDir, subject, expDate, num2str(expr), 'Episode001.h5');
lineNames = {'frameTrigger','FrameOut','CaptureActive'}; %17/7/20
[syncDataOut] = LoadSyncEpisodeFunction(fnam, lineNames);
%OR [FTOSignalValid, stimulusSignalValid]=filterFTObyStimulus(fnam) ??

%ind = strcmp({Timeline.hw.inputs.name}, 'neuralFrames');

switch TSFieldName
    case 'FrameOut'
        infoInXml = loadThorlabsExperimentXml(fullfile(TSDir, expt.subject, expt.expDate, ...
            num2str(expt.expNum)));
        
        if manualSave
            tgtData = filterFTObyManual(syncDataOut, infoInXml);%, 1);
        else
            [FTOSignalValid, stimulusSignalValid] = ...
                filterFTObyStimulusFastZ(syncDataOut, infoInXml, 1);% to show figures
            tgtData = FTOSignalValid;% .* stimulusSignalValid;
        end
        
        
        %thr = 0.35; % Set this lower as in some cases there are slow changes in phd amplitude. Check this if the number of frames seems incorrect
        thr = 0.5;
        
        above = tgtData>thr;
        deltas = [0; diff(above)];
        
        idx = find(deltas==1); %onset of the signal going HIGH
        
        nTotalFrames = length(idx);
        
    case 'FrameCounter'
        %if isfield(syncDataOut, 'FrameCounter') %not yet tested
        tgtData = syncDataOut.FrameCounter;
        
        nTotalFrames = tgtData(end);
        TTLs = [0; diff(tgtData)];
        idx = find(TTLs);
        
end

if isequal(frames, 'all')
    frames = 1:nTotalFrames;
else
    frames = frames(frames<=nTotalFrames);
end

times_all = syncDataOut.time(idx(frames));

TIdata_all = [];
if isfield(syncDataOut, 'PiezoMonitor')
    TIdata_all.PiezoMonitor = syncDataOut.PiezoMonitor(idx(frames));
end
if isfield(syncDataOut, 'PockelsMonitor')
    TIdata_all.PockelsMonitor = syncDataOut.PockelsMonitor(idx(frames));
end
if isfield(syncDataOut, 'FrameCounter')
    TIdata_all.FrameCounter = syncDataOut.FrameCounter(idx(frames));
end
if isfield(syncDataOut, 'npi')
    TIdata_all.npi = syncDataOut.npi(idx(frames));
end

%% for fastz acquisition mode

if infoInXml.numZSlices > 1
    movStride = infoInXml.numZSlices;% + infoInXml.numFlybackFrames; %10/2/20
else
    movStride = 1;
end
times = [];
TIdata = [];
for iplane = 1:infoInXml.numZSlices
    times(:,iplane) = times_all(iplane:movStride:end);
    
    if isfield(syncDataOut, 'PiezoMonitor')
        TIdata(:,iplane).PiezoMonitor = TIdata_all.PiezoMonitor(iplane:movStride:end);
    end
    if isfield(syncDataOut, 'PockelsMonitor')
        TIdata(:,iplane).PockelsMonitor = TIdata_all.PockelsMonitor(iplane:movStride:end);
    end
    if isfield(syncDataOut, 'FrameCounter')
        TIdata(:,iplane).FrameCounter = TIdata_all.FrameCounter(iplane:movStride:end);
    end
     if isfield(syncDataOut, 'npi')
        TIdata(:,iplane).npi = TIdata_all.npi(iplane:movStride:end);
    end
end




