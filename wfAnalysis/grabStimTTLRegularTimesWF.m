function expt = grabStimTTLRegularTimesWF(expt, timeline_filename, ISI,  verbose)
%times = grabStimTTLRegularTimesWF(expt, timeline_filename, ISI,  verbose)
%returns times of laser onset from Timeline data (allow multiple laser onsets per trial)
%note this function does not read protocol file, so the order of output
%variables are not in the order stored in p-file
%
% INPUT:
% expt.subject, expt.expDate, expt.expNum
% ISI: inter-stimulus interval[s]. used to detect stimulus onset
%
% OUTPUT:
% expt.DMDTimes.trialOn: start time of each trial
% expt.DMDTimes.onset: time when DMD image is switched AND laser was on
% expt.DMDTimes.laserPower: laser power of each image in [V]
%
% created from grabLaserTimesWF.m

% TODO:
% obviate expt.DMDTimes.trialOn(probably identical to grabStimTimesWF)
% check if the function works for multiple conditions

if nargin < 4 || isempty(verbose)
    verbose = 0; % Display progress or not
end

if nargin < 3 || isempty(ISI)
    ISI = 1/5; %[s] inter-repeat interval
end

if nargin < 2 || isempty(timeline_filename)
    thisDate = expt.expDate(1:10);
    thisSeries = str2num(expt.expDate(12:end));
    timeline_filename = dat.expFilePath(expt.subject, thisDate, thisSeries, expt.expNum, 'Timeline', 'master');
end


% filename_TL = sprintf('%s_%d_%s_Timeline.mat',expt.expDate, expt.expNum, expt.subject);
% timeline_filename = fullfile(TLDir, expt.subject, expt.expDate, num2str(expt.expNum), filename_TL);

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
    
    
    %% output from DMD
    DMD_idx = strcmp({Timeline.hw.inputs.name}, 'DMDOut');%'laserIn');%'DMD0');
    DMDTh = 2.5;%mean(Timeline.rawDAQData(:,DMD_idx)); %25/9/20
    DMD_trace_raw = Timeline.rawDAQData(:,DMD_idx)';
    DMD_on = DMD_trace_raw > DMDTh;
    
    DMD_trace = zeros(size(DMD_trace_raw));
    DMD_trace(DMD_on)=1; %[0 1]
    
    %this judges whether the DMD signal flipped or continued
    DMD_flip = find((~DMD_trace(1:end-1) & DMD_trace(2:end)) | ...
        (DMD_trace(1:end-1) & ~DMD_trace(2:end)))+1; %time idx when signal flips 0 to 1 or 1 to 0 ??
    
    
    DMD = struct('timestamps',[],'values',[]);
    DMD.timestamps = tltime(DMD_flip);
    DMD.values = DMD_trace(DMD_flip)';
    
    % DMD_offsets_time = DMD.timestamps(DMD.values == 0); %maybe missing/adding last frame of each repeat
    DMD_onsets_time = DMD.timestamps(DMD.values == 1)';
    DMD_offsets_time = DMD.timestamps(DMD.values == 0)';
    DMD_onset_tlidx = DMD_flip(DMD.values==1);
    DMD_offset_tlidx = DMD_flip(DMD.values==0);
    
    %% laser trace
    laserIn_idx = strcmp({Timeline.hw.inputs.name}, 'laserIn');
    %laserTh = 0;
    
    laser_raw = Timeline.rawDAQData(:,laserIn_idx)';
    %     laser_on = laser > lsrTh;
    %     %more robust?: laser_on = laser_on .* stimScreen_on;
    %     srate = 1/median(diff(tltime));
    %     laser_on_t = tltime(laser_on);
    
    % if isempty(lsrTh)
    %     laser_tmp = laser_raw;
    %     laser_tmp(DMD_trace==0) = 0;
    %     lsrTh = median(laser_tmp(laser_tmp>0));
    % end
    % 
    % laser = zeros(size(DMD_trace)); %[0 1];
    % laser(laser_raw > lsrTh) = 1;
    % laser(DMD_trace==0) = 0; 
    % laser_t = find(laser);
    
 
    ntrials = length(syncSquare_white_start_t);
    % laserOn_times = nan(ntrials,1);
    % laserOff_times = nan(ntrials,1);
    DMDOn_times = [];
    laserPower = [];
    for tt = 1:ntrials
        thisTrial = find((syncSquare_white_start_t(tt)- ISI/2 <= DMD_onsets_time) .* ...
            (syncSquare_white_end_t(tt) + ISI/2 >= DMD_onsets_time));
        % [DMD_off_idx] = find((syncSquare_white_start_t(tt)- ISI/2 <= DMD_offsets_time) .* ...
        %     (syncSquare_white_end_t(tt) + ISI/2 >= DMD_offsets_time));

        % remove events when DMD image is replaced but no laser output (at the start/end of each trial)
        lsrTh = 0.1;
        laserPower_thisTrial = [];
        for ii = 1:numel(thisTrial)
            laserPower_thisTrial(ii) = mean(laser_raw(DMD_onset_tlidx(thisTrial(ii)):DMD_offset_tlidx(thisTrial(ii))));
        end
        laserOut = find( laserPower_thisTrial > lsrTh);
        thisTrial = thisTrial(laserOut);

        if ~isempty(thisTrial)
            %laserOn_times(tt) = DMD_onsets(DMD_on_idx);  %not laser on
            %laserOff_times(tt) = DMD_offsets(DMD_off_idx);  %not laser on

            DMDOn_times = [DMDOn_times DMD_onsets_time(thisTrial)];
            %DMDOff_times = DMD_offsets(DMD_off_idx);

            laserPower = [laserPower laserPower_thisTrial(laserOut)];
        end
    end
    expt.DMDTimes.trialOn = syncSquare_white_start_t;
    expt.DMDTimes.onset = DMDOn_times;
    % expt.DMDTimes.offset = Off_times;
    expt.DMDTimes.laserPower = laserPower;
end
