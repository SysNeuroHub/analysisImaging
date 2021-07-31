function times = getStimTimesWF(expt, getAllPhdFlips, TLDir, IRI)
%times = getStimTimesWF
%returns times of visual stimulus onset from Timeline data
%
% INPUT:
% expt.subject, expt.expDate, expt.expNum
% getAllPhdFlips: if 1, return times.frameTimes
% IRI: inter-repeat interval[s]. used to detect repeat onset
%
% OUTPUT:
% times.onset
% times.offset
% times.frameTimes{iStim}
%
% assumes flicker SyncSquare & following timeline entries:
% photoDiode ... measured photodiode signal from screen
% syncSquare ... copy of syncsquare (without backlight ctrl) (from vs pc)
% stimScreen ... signal sent to backlight control (from arduino)
%
% 30/4/20 DS created from lilrig_load_experiment.m

% TODO: merge this with grabStimTimes for ts/ti

if nargin < 4
    IRI = 1/5; %[s] inter-repeat interval
end
verbose = false;


%% Display progress or not
if ~exist('verbose','var')
    verbose = false;
end

%% Load timeline

% [timeline_filename, timeline_exists] = AP_cortexlab_filename(expt.subject,expt.expDate,expt.expNum,'timeline');
% TLDir = fullfile('C:\Users\Experiment\Documents\MATLAB\Data\Subjects\',...
%     expt.subject,expt.expDate,num2str(expt.expNum));
filename_TL = sprintf('%s_%d_%s_Timeline.mat',expt.expDate, expt.expNum, expt.subject);
timeline_filename = fullfile(TLDir, filename_TL);
timeline_exists = 1;

if ~timeline_exists
    error([expt.subject ' ' expt.expDate ': no timeline']);
end

if timeline_exists
    if verbose; disp('Loading timeline...'); end
    
    load(timeline_filename);
    
    % Get wheter syncSquare on screen was flickering
    syncSquare_idx = strcmp({Timeline.hw.inputs.name}, 'syncSquare');
    if any(syncSquare_idx)
        syncSquare_flicker = max(Timeline.rawDAQData(:,syncSquare_idx)) - ...
            min(Timeline.rawDAQData(:,syncSquare_idx)) > 2;
        syncSquare_thresh = max(Timeline.rawDAQData(:,syncSquare_idx))/2;
        syncSquare_white_on = Timeline.rawDAQData(:,syncSquare_idx) > syncSquare_thresh; %[0 1]
        syncSquare_white_on_t = Timeline.rawDAQTimestamps(syncSquare_white_on);
        
        white_off = find(diff(syncSquare_white_on_t) > 2/55);% assume screen resolution at 60Hz
        syncSquare_white_start_t = [syncSquare_white_on_t(1) syncSquare_white_on_t(white_off + 1)]; % ~stimulus start
        syncSquare_white_end_t = [syncSquare_white_on_t(white_off) syncSquare_white_on_t(end)]; % ~stimulus end
        syncSquare_withinStim = logical(zeros(numel(Timeline.rawDAQTimestamps),1));
        for irepeat = 1:numel(syncSquare_white_start_t)
            syncSquare_withinStim((Timeline.rawDAQTimestamps>=syncSquare_white_start_t(irepeat)) & ...
                (Timeline.rawDAQTimestamps<=syncSquare_white_end_t(irepeat))) = 1;
        end
        syncSquare_withinStim = logical(syncSquare_withinStim);
        syncSquare_black_on = logical(syncSquare_withinStim .* ~syncSquare_white_on);
        syncSquare_black_on_t = Timeline.rawDAQTimestamps(syncSquare_black_on);
    end
    
    % Get whether backlight was flickering
    % this is equal to input to backlight control
    % this channel is not very useful??
    stimScreen_idx = strcmp({Timeline.hw.inputs.name}, 'stimScreen');
    if any(stimScreen_idx)
        stimScreen_flicker = max(Timeline.rawDAQData(:,stimScreen_idx)) - ...
            min(Timeline.rawDAQData(:,stimScreen_idx)) > 2;
        stimScreen_thresh = max(Timeline.rawDAQData(:,stimScreen_idx))/2;
        stimScreen_on = Timeline.rawDAQData(:,stimScreen_idx) > stimScreen_thresh; %[0 1]
    end
    
    
    
    % Get flicker or steady photodiode
    %     photodiode_type = myScreenInfo.SyncSquare.Type;
    % THE CODE BELOW SHOULD WORK IN BOTH STEADY AND FLICKER BUT NOT YET CHECKED
    
    photodiode_idx = strcmp({Timeline.hw.inputs.name}, 'photoDiode');
    %     switch lower(photodiode_type)
    %         case 'flicker' %assumes black background (RigInfo.BackgroundColor=[0 0 0])
    %when backlight is off. complete darkness
    bkOff = mode(Timeline.rawDAQData(~stimScreen_on,photodiode_idx));
    
    %when backlight is on and stimulus is NOT presented ...
    bkOn_nostim = mode(Timeline.rawDAQData(logical(~syncSquare_withinStim.*stimScreen_on),photodiode_idx));
    
    screenOnTh = ( bkOff + bkOn_nostim)/2;
    
    photodiode_on = Timeline.rawDAQData(:,photodiode_idx) > screenOnTh;
    %more robust?: photodiode_on = photodiode_on .* stimScreen_on;
    
    %time when photodiode detect backlight on
    photodiode_on_t = Timeline.rawDAQTimestamps(photodiode_on);
    %             photodiodeOn_filt = Timeline.rawDAQData(photodiode_on, photodiode_idx); %medfilt1(Timeline.rawDAQData(stimScreen_on, photodiode_idx),5);
    %             photodiode_thresh = median(photodiodeOn_filt);
    
    % 1 when flicker syncsquare is HIGH, 0 when flicker syncsquare is LOW
    %             photodiode_trace = photodiodeOn_filt > photodiode_thresh;
    photodiode_trace = syncSquare_white_on(photodiode_on); %4/5 use syncSquare w/b!
    
    %this judges whether the photodiode signal flipped or continued
    photodiode_flip = find((~photodiode_trace(1:end-1) & photodiode_trace(2:end)) | ...
        (photodiode_trace(1:end-1) & ~photodiode_trace(2:end)))+1; %time idx when signal flips 0 to 1 or 1 to 0 ??
    
    
    photodiode = struct('timestamps',[],'values',[]);
    photodiode.timestamps = photodiode_on_t(photodiode_flip)';
    photodiode.values = photodiode_trace(photodiode_flip);
    %< this is used to detect:
    % 1, stimulus on/offset of each repeat
    % 2, onset of each stimulus frame, specified in vs
    
    
    %% sanity check figures
    if verbose
        figure;
        ax(1) = subplot(311);
        plot(Timeline.rawDAQTimestamps, Timeline.rawDAQData(:,photodiode_idx), ...
            photodiode_on_t, Timeline.rawDAQData(photodiode_on, photodiode_idx));
        grid on;
        line(Timeline.rawDAQTimestamps([1 end]), [screenOnTh screenOnTh],'color','k')
        legend('raw PD', 'raw PD resample','th');
        
        ax(2) = subplot(312);
        plot(Timeline.rawDAQTimestamps,  Timeline.rawDAQData(:,syncSquare_idx));
        grid on;
        legend('syncSquare');
        
        ax(3) = subplot(313);
        plot(photodiode_on_t, photodiode_trace, '.-',...
            photodiode.timestamps, photodiode.values, '*'); %....looks like this is miraculuously correct
        grid on;
        legend('PD trace', 'PD flip');
        ylim([-0.2 1.2])
    end
    
    %         case 'steady'
    %             warning('dont use Steady syncSquare!!');
    %     end
    photodiode_offsets = photodiode.timestamps(photodiode.values == 0); %maybe missing/adding last frame of each repeat
    photodiode_onsets = photodiode.timestamps(photodiode.values == 1);
    
    %ANOTHER WAY TO DETECT ON/OFFSETS. IS THIS MORE ROBUST??
    % Stim times should just be odd (on) and even (off)
    %     if mod(length(photodiode_flip_times),2) == 0
    %         photodiode_onsets = photodiode_flip_times(1:2:end);
    %         photodiode_offsets = photodiode_flip_times(2:2:end);
    %     else
    %         error('Odd number of photodiode flips')
    %     end
    
    % Get specific stim onsets by time between last offset and new onset
    % (occasionally there a bad frame so flip but not new stim)
    stimOn_idx = [1;find(photodiode_onsets(2:end) - photodiode_offsets(1:end-1) > IRI) + 1];
    stimOff_idx = [stimOn_idx(2:end)-1; length(photodiode_offsets)];
    stimOn_times = photodiode_onsets(stimOn_idx);
    stimOff_times = photodiode_offsets(stimOff_idx);
    if verbose
        ax(3) = subplot(313);
        hold on;
        plot(stimOn_times, 0, 'o', stimOff_times, 0, 'x');
        linkaxes(ax,'x');
    end
end

times.onset = stimOn_times;
times.offset = stimOff_times; %maynot be accurate
if getAllPhdFlips
    for iStim = 1:length(stimOn_times)
        tmpUp = photodiode_onsets(photodiode_onsets>=times.onset(iStim) ...
            & photodiode_onsets<=times.offset(iStim));
        tmpDown = photodiode_offsets(photodiode_offsets>=times.onset(iStim) ...
            & photodiode_offsets<=times.offset(iStim));
        times.frameTimes{iStim} = sort([tmpUp(:); tmpDown(:)], 'ascend');
    end
end
