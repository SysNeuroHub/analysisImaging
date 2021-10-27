 function [Time, StartTime, EndTime] = ...
    alignTimeTL(Timeline, animal, iSeries, iExp, iStim, iTr, inputNameTL, pcoedge)
% Time = alignTimeTL(Timeline, animal, iSeries, iExp, iStim, iTr)
% returns time of signal on/offset relative to Timeline recording onset.
% If iStim = [] & iTr = [], then returns Time of all the time stamps
% recorded in Timeline
%
% Time = alignTimeTL(Timeline, animal, iSeries, iExp, iStim, iTr, inputNameTL)
% specifies name in timeline data
%
% Time = alignTimeTL(Timeline, animal, iSeries, iExp, iStim, iTr, inputNameTL, pcoedge)
% is for pco.edge camera only
%
% INPUT:
% Timeline: Timeline structure from mpep / choiceWorld experiment
% inputNameTL: Timeline.hw.inputs.name for camera input (e.g. cam1) (default: cam2)
% pcoedge (ScanMode, ShutterMode, FrameRate, nRows_raw)

% % currently not usable for choice world recording (only mpep experiment)
% 2014-03-05 DS created
% 2014-04-28 DS separated part of this file to extractTimeTL.m and extractSequenceTL.m
% 2014-08-05 DS added 8th input for pco.edge camera
% 2015-07-24 DS now time is from Timeline onset (not from the vis stim
% onset)
% 2016-02-12 DS allow iseries and iexp to be empty, in this case all the
% record in Timeline is used

% TO DO
% convert 2nd-4th input to Exps structure
% make compatible with flickering syncsquare?
% add stTime and enTime inputs, to be compatible with choice world data?
% remove pcoedge from input. obtain it inside this function
% remove animal, iSeries, iExp, iStim (and Timeline). replace it with
% timeline

doVisualize = false;

if nargin > 7 && ~isempty(pcoedge)
    if strcmpi(pcoedge.ScanMode, 'slow')
        trow = 27.52;
    elseif strcmpi(pcoedge.ScanMode, 'fast')
        trow = 9.17;
    else
        error('pcoedge.ScanMode was not recognized');
    end
end
 
if nargin < 7 || isempty(inputNameTL)
    inputNameTL = 'cam2';
end

if isempty(iStim) && isempty(iTr)
    stIdx = 1;
    enIdx = Timeline.rawDAQSampleCount;
else
    %find time index of stimStart/End
    [stIdx, enIdx] = tl.extractTimeTL(Timeline, animal, iSeries, iExp, iStim, iTr);
end
%vsStTime = Timeline.rawDAQTimestamps(vsStIdx);
%4/8/14 assuming 60Hz refresh rate

%extract data sequence (name specified by inputNameTL)
[dataSequence] = tl.extractSequenceTL(Timeline, stIdx, enIdx, inputNameTL);

%know time of exposure relative to visual stimulus onset
%time index for camera acquisition
if isequal(unique(diff(dataSequence)), [0 1]') %for the case of counter input
    Idx = find(diff(dataSequence) > 0) + stIdx;
    
else %for the case of analogue input
    aiTh = 0.5*max(dataSequence);
    PeriodIdx = find(dataSequence > aiTh) + stIdx;
    
    [EndIdx_cache] = find(diff(PeriodIdx) > 1);
    %EndIdx = EndIdx_cache + stIdx;
    EndIdx = [PeriodIdx(EndIdx_cache); max(PeriodIdx)];
    
    StartIdx = [min(PeriodIdx); PeriodIdx(EndIdx_cache + 1)];
    
    
    if length(EndIdx) < 10 %hack for M150819_SD_cam2
        aiTh = 0.75 * max(dataSequence);
        PeriodIdx = find(dataSequence > aiTh) + stIdx;
        
        [EndIdx_cache] = find(diff(PeriodIdx) > 1);
        %EndIdx = EndIdx_cache + stIdx;
        EndIdx = [PeriodIdx(EndIdx_cache); max(PeriodIdx)];
        
        StartIdx = [min(PeriodIdx); PeriodIdx(EndIdx_cache + 1)];
    end
    
    if nargin > 7 
        %disentangle merged frames. for M130518_SD
        jitterSize = 0.1; %steps in TL
        IdxDif = (EndIdx - StartIdx);
        typicalIdxDif = median(IdxDif);
        DifFromTypical = min(mod(IdxDif, typicalIdxDif), typicalIdxDif - mod(IdxDif, typicalIdxDif));
        mergedIdx = intersect(find(round(IdxDif/typicalIdxDif) >= 2), find(DifFromTypical < jitterSize*typicalIdxDif));
        %modified 6/5/14. IdxDif/typicalIdxDif should be equal or larger than 2
        
        
        if doVisualize
            plot(dataSequence);hold on;
            plot(PeriodIdx-stIdx, aiTh*ones(length(PeriodIdx),1),'w.');
            plot(StartIdx-stIdx, ones(length(StartIdx),1), 'ro'); % detected step onset
        end
        
        if ~isempty(mergedIdx) %TO DO: separated into 3 when two successive frames were merged?
            StartIdx = unique(sort([StartIdx; StartIdx(mergedIdx) + typicalIdxDif], 'ascend'));
            EndIdx = unique(sort([EndIdx; EndIdx(mergedIdx) - typicalIdxDif], 'ascend'));
            
            display([num2str(length(mergedIdx)) ' merged frames were separated (timeline).']);
            
            %to check
            if doVisualize
                plot(StartIdx(mergedIdx)-stIdx + typicalIdxDif, ...
                    ones(length(mergedIdx),1), 'm*') %disentangled step onset
                plot(EndIdx-stIdx, ones(length(EndIdx),1), 'bo'); % detected step offset
                plot(EndIdx(mergedIdx)-stIdx + typicalIdxDif, ...
                    ones(length(mergedIdx),1), 'c*') %disentangled step onset
            end
        end
    end
end

if isempty(StartIdx)
    error([ TLInputName ': Onset timing was not detected!']);
end
%need to check size of idx is consistent with mmapfile data


if nargin < 8 %% generic timeline input
    StartTime = Timeline.rawDAQTimestamps(StartIdx);
    %StartTime_fromvsSt = StartTime - vsStTime;
    EndTime = Timeline.rawDAQTimestamps(EndIdx);
    %EndTime_fromvsSt = EndTime - vsStTime;
   
    Time = 0.5*(StartTime + EndTime);
    %Time_fromvsSt = 0.5*(StartTime_fromvsSt + EndTime_fromvsSt);
    
else %% specific to pco.edge camera
    t_exposure = Timeline.rawDAQTimestamps(EndIdx+1) - Timeline.rawDAQTimestamps(StartIdx);
    
    if strcmpi(pcoedge.ShutterMode, 'rolling')
        t_frame = 0.5 * pcoedge.nRows_raw * trow/1e6; %frame readout time[s]
        
        %midest time of exposure at the middest pixel in rows
        Time = Timeline.rawDAQTimestamps(StartIdx) + t_exposure/2 + t_frame/2;
        
        if nargout >= 2
            StartTime = Timeline.rawDAQTimestamps(StartIdx);
            EndTime = Timeline.rawDAQTimestamps(EndIdx);
        end
        
        
        %% if timestamp is aggregated ex. M140522_LFR
        if prctile(t_exposure, 90) > 0.5
            display(['Each frame of ' inputNameTL ' was not detected. Use onset/offset to estimate timestamp.'])
            
            %rough approximation
            t_exposure = 1/pcoedge.FrameRate;
            
            %onset/offset time of the middle pixel in rows in Timeline [s]
            trOnTime = Timeline.rawDAQTimestamps(StartIdx(1)) + t_exposure/2 + t_frame/2;
            trOffTime = Timeline.rawDAQTimestamps(EndIdx(end)) + t_exposure/2 + t_frame/2;
            
            nFrames = round(pcoedge.FrameRate * ( trOffTime - trOnTime));%should add frate to input
            
            StartTime = (0:nFrames-1)/pcoedge.FrameRate + trOnTime;
            EndTime = (1:nFrames)/pcoedge.FrameRate + trOnTime;
            Time = 0.5*(StartTime + EndTime);
        end
        
    elseif strcmpi(pcoedge.ShutterMode, 'global')
        error('Global shutter mode is not yet implemented.');
    else
        error('Camera shutter mode not recognized.');
    end
end
    