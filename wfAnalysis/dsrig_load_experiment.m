
% 30/4/20 created from lilrig_load_experiment.m
% Loads data from experiments on dsrig

addpath('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\master\visbox\Stimulus'); %screenInfo

% Options to set before running:
%
% expt.subject - 'animal_name'
% expt.expDate - 'yyyy-mm-dd_s'
% expt.expNum - experiment number
% load_parts.cam/imaging/ephys = true/false (specify what to load,
% everything by default)
% verbose = true/false (display progress in command line, false default)

expt.subject = 'dummy_wf';
expt.expDate = '2020-04-30_1';
% expt.expNum = 2;%bklight ctr: ON,  syncSquare: steady,  camera 35Hz
% expt.expNum = 3;%bklight ctr: ON,  syncSquare: flicker, camera 35Hz
% expt.expNum = 4;%bklight ctr: OFF, syncSquare: flicker, camera 35Hz
% expt.expNum = 10;%bklight ctr: ON, syncSquare: steady, camera 70Hz
% expt.expNum = 11;%bklight ctr: ON, syncSquare: flicker, camera 70Hz
% expt.expNum = 14;%bklight ctr: ON, syncSquare: flicker, camera 70Hz, syncSquare recorded
 expt.expNum = 15;%bklight ctr: OFF, syncSquare: flicker, camera 70Hz, syncSquare recorded


%% Display progress or not
if ~exist('verbose','var')
    verbose = false;
end

%% Define what to load


% If nothing specified, load everything
if ~exist('load_parts','var')
    load_parts.cam = true;
    load_parts.imaging = true;
    load_parts.ephys = true;
else
    % If only some things specified, don't load others
    if ~isfield(load_parts,'cam')
        load_parts.cam = false;
    end
    if ~isfield(load_parts,'imaging')
        load_parts.imaging = false;
    end
    if ~isfield(load_parts,'ephys')
        load_parts.ephys = false;
    end
end

%% Load timeline and associated inputs

% [timeline_filename, timeline_exists] = AP_cortexlab_filename(expt.subject,expt.expDate,expt.expNum,'timeline');
loadDir_TL = fullfile('C:\Users\Experiment\Documents\MATLAB\Data\Subjects\',...
    expt.subject,expt.expDate,num2str(expt.expNum));
filename_TL = sprintf('%s_%d_%s_Timeline.mat',expt.expDate, expt.expNum, expt.subject);
timeline_filename = fullfile(loadDir_TL, filename_TL);
timeline_exists = 1;

% Load in protocol
% [protocol_filename, protocol_exists] = AP_cortexlab_filename(expt.subject,expt.expDate,expt.expNum,'protocol');
load(fullfile(loadDir_TL,'Protocol.mat')); %created in mpep pc

% Load in hardware
%HWexpPath = [fileparts(TLexpPath), sprintf('\\%s_%s_%s_hardwareInfo.mat', day, experiment, mpep_animal)];
filename_HW = sprintf('%s_%d_%s_hardwareInfo.mat',expt.expDate, expt.expNum, expt.subject);
load(fullfile(loadDir_TL, filename_HW));


if ~timeline_exists
    error([expt.subject ' ' expt.expDate ': no timeline']);
end

if timeline_exists
    if verbose; disp('Loading timeline...'); end
    
    load(timeline_filename);
    
    % Get camera times
    cam_name = 'camExposure';%'pcoExposure';
    timeline_cam_idx = strcmp({Timeline.hw.inputs.name}, cam_name);
    
    cam_expose_starts = Timeline.rawDAQTimestamps( ...
        find(Timeline.rawDAQData(1:end-1,timeline_cam_idx) <= 2 & ...
        Timeline.rawDAQData(2:end,timeline_cam_idx) > 2) + 1);
    cam_expose_stops = Timeline.rawDAQTimestamps( ...
        find(Timeline.rawDAQData(1:end-1,timeline_cam_idx) >= 2 & ...
        Timeline.rawDAQData(2:end,timeline_cam_idx) < 2) + 1);
    
    cam_time = cam_expose_starts;
    cam_expose_times = cam_expose_stops - cam_expose_starts;
    
    % Get acqLive signal
    acqLive_name = 'acqLive';
    acqLive_idx = strcmp({Timeline.hw.inputs.name}, acqLive_name);
    thresh = max(Timeline.rawDAQData(:,acqLive_idx))/2;
    acqLive_trace = Timeline.rawDAQData(:,acqLive_idx) > thresh;
    acqLive_timeline = Timeline.rawDAQTimestamps( ...
        [find(acqLive_trace,1),find(acqLive_trace,1,'last')+1]);
    
    %     % Get flipper signal (this was added late, might not be present)
    %     % AP records the same in one ephys channel, later match ephys and timeline
    %     flipper_name = 'flipper';
    %     flipper_idx = strcmp({Timeline.hw.inputs.name}, flipper_name);
    %     flipper_thresh = 2; % TTL threshold
    %     flipper_trace = Timeline.rawDAQData(:,flipper_idx) > flipper_thresh;
    %     flipper_flip = find((~flipper_trace(1:end-1) & flipper_trace(2:end)) | ...
    %         (flipper_trace(1:end-1) & ~flipper_trace(2:end)))+1;
    %     flipper_flip_times_timeline = Timeline.rawDAQTimestamps(flipper_flip)';
    
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
    photodiode_type = myScreenInfo.SyncSquare.Type;
    
    photodiode_idx = strcmp({Timeline.hw.inputs.name}, 'photoDiode');
    switch lower(photodiode_type)
        case 'flicker' %assumes black background (RigInfo.BackgroundColor=[0 0 0])
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
                linkaxes(ax,'x');
            end
            
        case 'steady'
            warning('dont use Steady syncSquare!!');
    end
    photodiode_offsets = photodiode.timestamps(photodiode.values == 0); %maybe missing/adding last frame of each repeat
    photodiode_onsets = photodiode.timestamps(photodiode.values == 1);
    
    %IS THE BELOW CORRECT?
     % Stim times should just be odd (on) and even (off)
     %     if mod(length(photodiode_flip_times),2) == 0
     %         photodiode_onsets = photodiode_flip_times(1:2:end);
     %         photodiode_offsets = photodiode_flip_times(2:2:end);
     %     else
     %         error('Odd number of photodiode flips')
     %     end
    
    % Get specific stim onsets by time between last offset and new onset
    % (occasionally there a bad frame so flip but not new stim)
    refresh_rate_cutoff = 1/5; %[s] inter-repeat interval
    stimOn_times = photodiode_onsets( ...
        [1;find(photodiode_onsets(2:end) - photodiode_offsets(1:end-1) > refresh_rate_cutoff) + 1]);
        
end



%% Load mpep protocol

% [protocol_filename,protocol_exists] = AP_cortexlab_filename(animal,day,experiment,'protocol');
% 
% if protocol_exists
    
%     if verbose; disp('Loading mpep protocol...'); end
%     
%     load(protocol_filename);
%     
%     % Load in hardware info
%     hwinfo_filename = AP_cortexlab_filename(animal,day,experiment,'hardware');
%     load(hwinfo_filename);
    
%     % Stim times should just be odd (on) and even (off)
%     if mod(length(photodiode_flip_times),2) == 0
%         photodiode_onsets = photodiode_flip_times(1:2:end);
%         photodiode_offsets = photodiode_flip_times(2:2:end);
%     else
%         error('Odd number of photodiode flips')
%     end
    
    % Get flicker/steady photodiode mode
    %photodiode_type = lower(myScreenInfo.SyncSquare.Type);
    
    % Get stim on times
%     if strcmp(Protocol.xfile,'stimSparseNoiseUncorrAsync.x')
%         % Sparse noise
%         % (at the moment this is in a sparse noise-specific script)
%     else
        %         % Anything else
        %         % Get specific stim onsets by time between last offset and new onset
        %         % (occasionally there a bad frame so flip but not new stim)
        %         refresh_rate_cutoff = 1/5;
        %         stimOn_times = photodiode_onsets( ...
        %             [1;find(photodiode_onsets(2:end) - photodiode_offsets(1:end-1) > refresh_rate_cutoff) + 1]);
        %
        %         if length(stimOn_times) ~= numel(Protocol.seqnums)
        %             error('MPEP/Photodiode error: photodiode doesn''t match stim')
        %         end
        
        stimIDs = zeros(size(stimOn_times));
        for q = 1:size(Protocol.seqnums,1)
            stimIDs(Protocol.seqnums(q,:)) = q;
        end
%     end
    
% end
