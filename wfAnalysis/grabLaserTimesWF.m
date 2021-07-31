function expt = grabLaserTimesWF(expt, TLDir, ISI,  verbose, AOTh, lsrTh)
%times = getStimTimesWF
%returns times of visual stimulus onset from Timeline data.
%note this function does not read protocol file, so the order of output
%variables are not in the order stored in p-file
%
% INPUT:
% expt.subject, expt.expDate, expt.expNum
% ISI: inter-stimulus interval[s]. used to detect stimulus onset
% analogSwitchCtrl: whether the laser is controled by arduino or not
% (default: true).
%
% OUTPUT:
% expt.stimTimes.onset: start time of each trial
%
% created from grabStimTimesWF.m
%
% TODO
% use "laser" stored in Ai7 for more precise detection of onset
% add:
% expt.stimTimes.offset: stop time of each trial
% expt.stimTimes.frameTimes{iStim}

if nargin < 6
    lsrTh = [];
end
if nargin < 4 || isempty(verbose)
    verbose = 0; % Display progress or not
end

if nargin < 3 || isempty(ISI)
    ISI = 1/5; %[s] inter-repeat interval
end

if nargin < 2 || isempty(TLDir)
    %thisDate = expt.expDate(1:10);
    %thisSeries = str2num(expt.expDate(12:end));
    %TLDir = fileparts(dat.expFilePath(expt.subject, thisDate, thisSeries, expt.expNum, 'timeline','master'));
    TLDir = '\\storage.erc.monash.edu.au\shares\MNHS-dshi0006\Subjects';
end


filename_TL = sprintf('%s_%d_%s_Timeline.mat',expt.expDate, expt.expNum, expt.subject);
timeline_filename = fullfile(TLDir, expt.subject, expt.expDate, num2str(expt.expNum), filename_TL);

if ~exist(timeline_filename,'file')
    error([expt.subject ' ' expt.expDate ': no timeline']);
else
    if verbose; disp('Loading timeline...'); end
    
    load(timeline_filename);
    
    tltime = Timeline.rawDAQTimestamps';
    
    %% Get trial onset times
    syncSquare_idx = strcmp({Timeline.hw.inputs.name}, 'syncSquare');
    if any(syncSquare_idx)
        syncSquare_flicker = max(Timeline.rawDAQData(:,syncSquare_idx)) - ...
            min(Timeline.rawDAQData(:,syncSquare_idx)) > 2;
        syncSquare_thresh = max(Timeline.rawDAQData(:,syncSquare_idx))/2;
        syncSquare_white_on = Timeline.rawDAQData(:,syncSquare_idx) > syncSquare_thresh; %[0 1]
        syncSquare_white_on_t = tltime(syncSquare_white_on);
        
        white_off = find(diff(syncSquare_white_on_t) > 2/55);% assume screen resolution at 60Hz
        syncSquare_white_start_t = [syncSquare_white_on_t(1); syncSquare_white_on_t(white_off + 1)]; % ~stimulus start
        syncSquare_white_end_t = [syncSquare_white_on_t(white_off); syncSquare_white_on_t(end)]; % ~stimulus end
        syncSquare_withinStim = logical(zeros(numel(tltime),1));
        for irepeat = 1:numel(syncSquare_white_start_t)
            syncSquare_withinStim((tltime>=syncSquare_white_start_t(irepeat)) & ...
                (tltime<=syncSquare_white_end_t(irepeat))) = 1;
        end
        syncSquare_withinStim = logical(syncSquare_withinStim);
        syncSquare_black_on = logical(syncSquare_withinStim .* ~syncSquare_white_on);
        syncSquare_black_on_t = tltime(syncSquare_black_on);
    end
    
    
    %% Get whether backlight was flickering
    % this is equal to input to backlight control
    % this channel is not very useful??
    stimScreen_idx = strcmp({Timeline.hw.inputs.name}, 'stimScreen');
    if any(stimScreen_idx)
        stimScreen_flicker = max(Timeline.rawDAQData(:,stimScreen_idx)) - ...
            min(Timeline.rawDAQData(:,stimScreen_idx)) > 2;
        %stimScreen_thresh = max(Timeline.rawDAQData(:,stimScreen_idx))/2;
        stimScreen_thresh = 2.5;%this way compatible with backlight is always on
        stimScreen_on = Timeline.rawDAQData(:,stimScreen_idx) > stimScreen_thresh; %[0 1]
    else
        stimScreen_on = ones(Timeline.rawDAQSampleCount,1);
    end
    
    
    %% AO trace
    ao_idx = strcmp({Timeline.hw.inputs.name}, 'ao0');
    %AOTh = mean(Timeline.rawDAQData(:,ao_idx)); %25/9/20
    ao_trace_raw = Timeline.rawDAQData(:,ao_idx)';
    ao_on = ao_trace_raw > AOTh;
    
    ao_trace = zeros(size(ao_trace_raw));
    ao_trace(ao_on)=1; %[0 1]
    
    %this judges whether the ao signal flipped or continued
    ao_flip = find((~ao_trace(1:end-1) & ao_trace(2:end)) | ...
        (ao_trace(1:end-1) & ~ao_trace(2:end)))+1; %time idx when signal flips 0 to 1 or 1 to 0 ??
    
    
    ao = struct('timestamps',[],'values',[]);
    ao.timestamps = tltime(ao_flip);
    ao.values = ao_trace(ao_flip)';
    
    ao_offsets = ao.timestamps(ao.values == 0); %maybe missing/adding last frame of each repeat
    ao_onsets = ao.timestamps(ao.values == 1);
    
    
    %% laser trace
    photodiode_idx = strcmp({Timeline.hw.inputs.name}, 'laser');
    %laserTh = 0;
    
    laser_raw = Timeline.rawDAQData(:,photodiode_idx)';
    %     laser_on = laser > lsrTh;
    %     %more robust?: laser_on = laser_on .* stimScreen_on;
    %     srate = 1/median(diff(tltime));
    %     laser_on_t = tltime(laser_on);
    
    if isempty(lsrTh)
        laser_tmp = laser_raw;
        laser_tmp(ao_trace==0) = 0;
        lsrTh = median(laser_tmp(laser_tmp>0));
    end
        
    laser = zeros(size(ao_trace)); %[0 1];
    laser(laser_raw > lsrTh) = 1;
    laser(ao_trace==0) = 0; 
    laser_t = find(laser);
    
    %     %laser_trace = syncSquare_white_on_adj(laser_on); %7/7/20
    %     %laser_trace = laser(laser_on); %NG
    %
    %
    %     %this judges whether the laser signal flipped or continued
    %         laser_flip = find((~laser_trace(1:end-1) & laser_trace(2:end)) | ...
    %             (laser_trace(1:end-1) & ~laser_trace(2:end)))+1; %time idx when signal flips 0 to 1 or 1 to 0 ??
    %     % < this does not work as the laser trace is noisy
    %
    %     laser = struct('timestamps',[],'values',[]);
    %     laser.timestamps = laser_on_t(laser_flip);
    %     laser.values = laser_trace(laser_flip)';
    %
    %     laser_offsets = laser.timestamps(laser.values == 0); %maybe missing/adding last frame of each repeat
    %     laser_onsets = laser.timestamps(laser.values == 1);
    %
    %     %heuristics numbers of onsets and offsets do not match 14/10/20
    %     if (numel(laser_onsets) == numel(laser_offsets) + 1)
    %         upDuration = laser_offsets(1:end)-laser_onsets(1:end-1);
    %         if isempty(find(upDuration < 0, 1))
    %             laser_offsets = [laser_offsets; laser_onsets(end)+median(upDuration)];
    %         end
    %     end
    %
    %     stimOn_idx = [1; find(diff(laser_onsets)>ISI) + 1]; %15/1/20 more robust than the above if long (200s) trial
    %
    %     stimOff_idx = [stimOn_idx(2:end)-1; length(laser_offsets)];
    %     stimOn_times = laser_onsets(stimOn_idx);
    %     stimOff_times = laser_offsets(stimOff_idx);
    
    ntrials = length(syncSquare_white_start_t);
    laserOn_times = nan(ntrials,1);
    laserOff_times = nan(ntrials,1);
    for tt = 1:ntrials
        [AO_on_idx] = find((syncSquare_white_start_t(tt)- ISI/2 <= ao_onsets) .* ...
            (syncSquare_white_end_t(tt) + ISI/2 >= ao_onsets));
        [AO_off_idx] = find((syncSquare_white_start_t(tt)- ISI/2 <= ao_offsets) .* ...
            (syncSquare_white_end_t(tt) + ISI/2 >= ao_offsets));
        if ~isempty(AO_on_idx)
            %laserOn_times(tt) = ao_onsets(AO_on_idx);  %not laser on
            %laserOff_times(tt) = ao_offsets(AO_off_idx);  %not laser on
            
            AOOn_times = ao_onsets(AO_on_idx);
            AOOff_times = ao_offsets(AO_off_idx);
            
            laserOnIdx = min(find(tltime(laser_t) >= AOOn_times));
            laserOffIdx = max(find(tltime(laser_t) <= AOOff_times));
            laserOn_times(tt) = tltime(laser_t(laserOnIdx)); %2/12/20
            laserOff_times(tt) = tltime(laser_t(laserOffIdx)); %2/12/20
            
        end
    end
    expt.laserTimes.trialOn = syncSquare_white_start_t;
    expt.laserTimes.onset = laserOn_times;
    expt.laserTimes.offset = laserOff_times;
end
