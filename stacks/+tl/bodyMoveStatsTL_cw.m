function [speed_median, speed_std, SPEED, TIME] = ...
    bodyMoveStatsTL_cw(Exps, stimList, repeatList, windowSizeSec, drawSummary)
%
% [speed_median, speed_std] = bodyMoveStatsTL(animal, iSeries, iExp, iStim)
% returns median and std of speed [m/s] (in a trial) of locomotion of every trials
% stimList: 1D array
% repeatList: 1D cell array
%
% [...] = bodyMoveStatsTL(Exps, iStim, repeatList)
% returns speed of specified trials
%
% [...] = bodyMoveStatsTL(Exps, iStim, repeatList, TLstring)
% allows to specify which Timeline channel to use (default:'rotaryEncoder')
%
%  [...] = bodyMoveStatsTL(Exps, iStim, repeatList, TLstring, windowSize)
% allows to pecify length of smoothing window in seconds
%
% [speed_median, speed_std, SPEED, TIME] = bodyMoveStatsTL(...) also returns cell
% array of speed of each trial, and time from vs onset

%SetDefaultDirs;

% 2014-10-21 DS : allow iStim to be a vector (stimlist)
% 2014-11-05 DS : allow repeatList to be a cell, specifying repeatList for
% each stimulus
% 2015-11-15 DS : now compatible with LEGO wheel used in choice world
% 2015-11-17 DS : inputs reorganized


if nargin < 4
    windowSizeSec = 0.05;
end

%p = ProtocolLoad(Exps);

%4/9/2015 DS for parallel computing
% p = ProtocolLoad(Exps.animal, Exps.iseries, Exps.iexp, 'donotload',...
%         '\\zserver.ioo.ucl.ac.uk\Data\trodes');

%15/11/2015
[p, block] = ProtocolLoadDS(Exps);

ppr = 100;%number of pulses per revolution %KÜBLER - 05.2400.1122.0100
ballDiameter = 0.062; %[m] diamter of LEGO wheel
mpp = pi*ballDiameter / ppr / 4; %[m/pulse] on treadmil %2014/11/28 fixed

if isempty(stimList)
    stimList = p.nstim;
end

if ~exist('repeatList','var')
    repeatList = [];
end
if isempty(repeatList) || nargin < 3
    repeatList = 1:p.nrepeats;%TO DO: change this for each exp
end

if ~iscell(repeatList)
    for jjj = 1:length(stimList)
        repeatList_cache{jjj} = repeatList;
    end
    repeatList = repeatList_cache;
end

%% load Timeline
seriesStr = num2str(Exps.iseries);
serieName = [seriesStr(1:4) '-' seriesStr(5:6) '-' seriesStr(7:8)];
tlname = sprintf('//zserver.ioo.ucl.ac.uk/Data/expInfo/%s/%s/%d/%s_%d_%s_Timeline.mat',...
    Exps.animal, ...
    serieName, ...
    Exps.iexp, ...
    serieName, Exps.iexp, Exps.animal);

disp('bodyMoveStatsTL:Loading Timeline...')
load(tlname);


%% load position of wheel
% in CW, wheel position is saved in block, but not in timeline
rotarySequence_all = block.inputSensorPositions;
time_all = block.inputSensorPositionTimes;
srate = 1./median(diff(time_all)); %[Hz]

%% 
[pcoedge.ShutterMode, pcoedge.ScanMode, pcoedge.FrameRate, pcoedge.nRows_raw] = ...
    tools.pcoedgeExtraInfo(Exps);


for jjj = 1:length(stimList)
    
    iStim = stimList(jjj);
    
    for iii = 1:length(repeatList{jjj})
        
        iTr = repeatList{jjj}(iii);
        
        [stIdx, enIdx] = ...
            tl.extractTimeTL(Timeline, Exps.animal, Exps.iseries, Exps.iexp, iStim, iTr);
        
        [~, stIdx_rot] = min(abs(time_all - Timeline.rawDAQTimestamps(stIdx)));
        [~, enIdx_rot] = min(abs(time_all - Timeline.rawDAQTimestamps(enIdx)));
        
        
        rotarySequence = rotarySequence_all(stIdx_rot:enIdx_rot);
        time_rot = time_all(stIdx_rot:enIdx_rot-1) + 0.5/srate;
        
        % this cannot be used as acquisition onset in choice world
%         vsStTime = tl.visStimStart(Timeline, Exps.animal, Exps.iseries, Exps.iexp, ...
%             iStim, iTr); %14/8/15


%% get acquisition start time
    timeVecTL = ...
        tl.alignTimeTL(Timeline, Exps.animal, ...
        Exps.iseries, Exps.iexp, iStim, iTr, 'cam2', pcoedge);
     
    timeInterval = diff(timeVecTL);%interval between successive camera records
    typicalInterval = median(timeInterval);
    acqStartIdx_cache = find(timeInterval > 10*typicalInterval) + 1;
    acqStartTime_cache = timeVecTL(acqStartIdx_cache);%onset(s) of camera acquisition within the extracted period
    
    if ~isempty(acqStartIdx_cache)
        [~, minidx] = min(abs(acqStartTime_cache - block.trial(iTr).trialStartedTime));
        acqStartTime = acqStartTime_cache(minidx);
        %acqStartIdx = find(timeVecTL == acqStartTime);
    else
        acqStartTime = timeVecTL(1);
        %acqStartIdx = 1;
    end
    
%     [StackDir, ThisFileName] = ...
%          tools.getDirectory( ServerDir, p, Exps.ResizeFac, 1, iTr, Cam, suffix); 
%      load(fullfile(StackDir, ['info' ThisFileName]), 'info');
% 
%     if info.timeInfo.elim(1)
%         acqStartTime = TimeVecTL(acqStartIdx+1);
%     end


    %%
        
        windowSize = round(windowSizeSec * srate);
        
        rotaryFilt = filtfilt(ones(1,windowSize)/windowSize,1,double(rotarySequence));
        speed = srate*diff(rotaryFilt); %[pulses/sec]
        speed = mpp * speed; %[m/pulse]*[pulses/sec]
        
        speed_std(iii) = std(speed);
        speed_median(iii) = median(speed);
        
        if nargout > 2
            SPEED{jjj, iii} = speed;
        end
        if nargout > 3
            TIME{jjj,iii} = time_rot - acqStartTime;
        end
    end
end