function expt = grabLaserTimesWF(expt, TLDir, lsrTh,  verbose)
%expt = grabLaserTimesWF
%returns times of laser onset from Timeline data (assumed to happen once
%per trial)
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
% expt.laserTimes.onset: first time when laser power exceed threshold each trial
% expt.laserTimes.offset: last time when laser power exceed threshold each trial
%
% created from grabStimTimesWF.m


if nargin < 4 || isempty(verbose)
    verbose = 0; % Display progress or not
end



if nargin < 2 || isempty(TLDir)
    if ispc
        TLDir = '\\storage.erc.monash.edu.au\shares\MNHS-dshi0006\Subjects';
    elseif isunix
        TLDir = '//mnt/dshi0006_market/Subjects';
    end
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


    %% laser trace
    photodiode_idx = strcmp({Timeline.hw.inputs.name}, 'laserIn');
    %laserTh = 0;

    laser_raw = Timeline.rawDAQData(:,photodiode_idx)'; %[V]
  
    if nargin <3 || isempty(lsrTh)
        laser_tmp = laser_raw;
        %laser_tmp(ao_trace==0) = 0;
        lsrTh = median(laser_tmp(laser_tmp>0));
    end

    laser = zeros(size(laser_raw)); %[0 1];
    laser(laser_raw > lsrTh) = 1;
    %laser(ao_trace==0) = 0;
    laser_tidx = find(laser);

    ntrials = length(syncSquare_white_start_t);
    laserOn_times = nan(ntrials,1);
    laserOff_times = nan(ntrials,1);
    for tt = 1:ntrials
        idx_thisTrial = find((tltime(laser_tidx)>= syncSquare_white_start_t(tt)).*(tltime(laser_tidx)<= syncSquare_white_end_t(tt)));
        if ~isempty(idx_thisTrial)
            laserOn_times(tt) = min(tltime(laser_tidx(idx_thisTrial)));
            laserOff_times(tt) = max(tltime(laser_tidx(idx_thisTrial)));
        end
    end
    expt.laserTimes.trialOn = syncSquare_white_start_t;
    expt.laserTimes.onset = laserOn_times;
    expt.laserTimes.offset = laserOff_times;
end
