% Sparse noise retinotopy
%
% This script loads in data from sparse noise experiment, then finds
% retinotopy using a bootstrapped method that looks empirically the best

% addpath('C:\Users\Experiment\npy-matlab');
% addpath(genpath('C:\Users\Experiment\widefield_ns'));
% addpath('C:\Users\Experiment\Documents\MATLAB\widefield\generalUtils'); %lpFilt_DS
% addpath('\\zserver\Code\stacks\');
% addpath(genpath('C:\Users\Experiment\Documents\MATLAB\dsbox'));
addpath(genpath('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\dsbox'));
addpath('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\dsInUCL\MATLAB\dsbox');
addpath(genpath('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\master\visbox\rigbox'));
addpath('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\master\visbox\Stimulus'); %screenInfo

analyseImage = 0;

% %% Define options
% addpath('\\zserver.cortexlab.net\Code\Stimulus')
% addpath('\\zserver.cortexlab.net\Data\xfiles')
% % Animal and experiment details
% mpep_animal = 'CR_DS1';
% day = '2018-11-22';
% experiment = '1';%
rig = 'lilrig'; % kilotrode or bigrig
cam_color_n = 2;
cam_color_signal = 'blue';
cam_color_hemo = 'purple';


%% Load experiment info

% Set rig-specific timeline names
switch rig
    case 'bigrig'
        cam_name = 'cam2';
        acqLive_name = 'acqLiveEcho';
    case 'lilrig'
        cam_name = 'pcoExposure';
        acqLive_name = 'acqLive';
end

% expRef = dat.constructExpRef(mpep_animal, day, str2double(experiment));
% TLexpPath = dat.expFilePath(expRef, 'Timeline', 'master');
TLexpPath = 'F:\UCL\MonashExpData\Data\CR_DS3\2018-12-07\1\2018-12-07_1_CR_DS3_Timeline.mat';
load(TLexpPath, 'Timeline')

% Load timeline
%timeline_filename = get_cortexlab_filename(mpep_animal,day,experiment,'timeline');
%load(timeline_filename);
timeline_sample_rate = Timeline.hw.daqSampleRate;

% Get camera times
timeline_cam_idx = strcmp({Timeline.hw.inputs.name}, cam_name);
cam_samples = find(Timeline.rawDAQData(1:end-1,timeline_cam_idx) <= 2 & ...
    Timeline.rawDAQData(2:end,timeline_cam_idx) > 2);
cam_time = cam_samples./timeline_sample_rate;

% Load in hardware
%HWexpPath = [fileparts(TLexpPath), sprintf('\\%s_%s_%s_hardwareInfo.mat', day, experiment, mpep_animal)];
HWexpPath = 'F:\UCL\MonashExpData\Data\CR_DS3\2018-12-07\1\2018-12-07_1_CR_DS3_hardwareInfo.mat';
load(HWexpPath)
%hwinfo_filename = get_cortexlab_filename(mpep_animal,day,experiment,'hardware');
%load(hwinfo_filename);

% Load in protocol
%protocolexpPath = [fileparts(TLexpPath), '\Protocol.mat'];
protocolexpPath = 'F:\UCL\MonashExpData\Data\CR_DS3\2018-12-07\1\Protocol.mat';
load(protocolexpPath)
%try
%    protocol_filename = get_cortexlab_filename(mpep_animal,day,experiment,'protocol','8digit');
%    load(protocol_filename);
%catch me
%    protocol_filename = get_cortexlab_filename(mpep_animal,day,experiment,'protocol','day_dash');
%    load(protocol_filename);
%end

% Get acqLive signal
acqLive_idx = strcmp({Timeline.hw.inputs.name}, acqLive_name);
thresh = max(Timeline.rawDAQData(:,acqLive_idx))/2;
acqLive_trace = Timeline.rawDAQData(:,acqLive_idx) > thresh;
acqLive_timeline = Timeline.rawDAQTimestamps( ...
    [find(acqLive_trace,1),find(acqLive_trace,1,'last')+1]);

% Get stimulus onsets and parameters
photodiode_idx = strcmp({Timeline.hw.inputs.name}, 'photoDiode');
thresh = max(Timeline.rawDAQData(:,photodiode_idx))/2;
photodiode_trace = Timeline.rawDAQData(:,photodiode_idx) > thresh;
photodiode_flip_c = find((~photodiode_trace(1:end-1) & photodiode_trace(2:end)) | ...
    (photodiode_trace(1:end-1) & ~photodiode_trace(2:end)))+1;

photodiode_c = struct('timestamps',[],'values',[]); %this is used only for detecting stimulus onset
photodiode_c.timestamps = Timeline.rawDAQTimestamps(photodiode_flip_c)';
photodiode_c.values = photodiode_trace(photodiode_flip_c);

photodiode_onsets = photodiode_c.timestamps(photodiode_c.values == 1);

refresh_rate_cutoff = 1/5;
stim_onsets = photodiode_onsets( ...
    [1;find(diff(photodiode_onsets) > refresh_rate_cutoff) + 1]);

stimIDs = zeros(size(stim_onsets));
for q = 1:size(Protocol.seqnums,1)
    stimIDs(Protocol.seqnums(q,:)) = q;
end

%%
% Get flicker or steady photodiode
photodiode_type = myScreenInfo.SyncSquare.Type;

% Get stimulus onsets and parameters
photodiode_idx = strcmp({Timeline.hw.inputs.name}, 'photoDiode');

% Get stim screen signal index
stimScreen_idx = strcmp({Timeline.hw.inputs.name}, 'stimScreen');
if any(stimScreen_idx)
    stimScreen_flicker = max(Timeline.rawDAQData(:,stimScreen_idx)) - ...
        min(Timeline.rawDAQData(:,stimScreen_idx)) > 2;
    stimScreen_thresh = max(Timeline.rawDAQData(:,stimScreen_idx))/2;
    stimScreen_on = Timeline.rawDAQData(:,stimScreen_idx) > stimScreen_thresh; %[0 1]
end

switch lower(photodiode_type)
    case 'flicker'
        warning('if flickering photodiode and steady screen, write diff')
        %         % get differential of photodiode
        %         photodiode_diff = [diff(Timeline.rawDAQData(:,photodiode_idx));0];
        %
        %         stimScreen_idx = strcmp({Timeline.hw.inputs.name}, 'stimScreen');
        %         if any(Timeline.rawDAQData(:,stimScreen_idx) < 1);
        %             stimScreen_thresh = max(Timeline.rawDAQData(:,stimScreen_idx))/2;
        %             stimScreen_on = Timeline.rawDAQData(:,stimScreen_idx) > stimScreen_thresh;
        %             photodiode_diff(~stimScreen_on) = NaN;
        %             photodiode_
        %
        %         end
        
        % This is the same as below... can probably just use
        %             stimScreen_on = Timeline.rawDAQData(:,photodiode_idx) > 0.15; %[0 1]
        %             stimScreen_on_t = Timeline.rawDAQTimestamps(stimScreen_on);
        
        %             photodiode_thresh = max(Timeline.rawDAQData(:,photodiode_idx))/2;
        %             % median filter because of weird effect where
        %             % photodiode dims instead of off for one sample
        %             % while backlight is turning off
        %             photodiode_trace = medfilt1(Timeline.rawDAQData(stimScreen_on, ...
        %                 photodiode_idx),5) > photodiode_thresh; %1: white rectangle, 0: black rectangle?
        
        
        %% DS ... looks like this part is very fragile!!
        
%         stimScreen_on = Timeline.rawDAQData(:,photodiode_idx) > 0.5; %0.15 %NEED TWEAK
%         stimScreen_on_t = Timeline.rawDAQTimestamps(stimScreen_on); %time when photodiode detect stimulus on screen
         bkOn_stim0 = mode(Timeline.rawDAQData(stimScreen_on,photodiode_idx)); 
            
             %screenOnTh = ( bkOff + bkOn_stim0)/2;
            screenOnTh = bkOn_stim0;
            stimScreen_on = Timeline.rawDAQData(:,photodiode_idx) > screenOnTh; %0.15
            stimScreen_on_t = Timeline.rawDAQTimestamps(stimScreen_on); %time when photodiode detect stimulus on screen
 
        photodiodeOn_filt = Timeline.rawDAQData(stimScreen_on, photodiode_idx); %medfilt1(Timeline.rawDAQData(stimScreen_on, photodiode_idx),5);
        photodiode_thresh = median(photodiodeOn_filt);
        photodiode_trace = photodiodeOn_filt > photodiode_thresh; % 1 when flicker syncsquare is HIGH, 0 when flicker syncsquare is LOW
        
        %  plot(Timeline.rawDAQTimestamps, Timeline.rawDAQData(:,photodiode_idx), ...
        %      stimScreen_on_t, Timeline.rawDAQData(stimScreen_on, photodiode_idx),...
        %      stimScreen_on_t, photodiodeOn_filt)
        %%%%
        
        %this judges whether the photodiode signal flipped or continued
        photodiode_flip = find((~photodiode_trace(1:end-1) & photodiode_trace(2:end)) | ...
            (photodiode_trace(1:end-1) & ~photodiode_trace(2:end)))+1; %time idx when signal flips 0 to 1 or 1 to 0 ??
        
        photodiode = struct('timestamps',[],'values',[]);
        photodiode.timestamps = stimScreen_on_t(photodiode_flip)';
        photodiode.values = photodiode_trace(photodiode_flip);
        
    case 'steady'
        
        % Take into account if the screen flickers
        
        % have to redefine periods of screen on, because
        % sometimes there's a sample or so difference
        stimScreen_on = Timeline.rawDAQData(:,photodiode_idx) > 0.3;
        stimScreen_on_t = Timeline.rawDAQTimestamps(stimScreen_on);
        photodiode_thresh = (max(Timeline.rawDAQData(:,photodiode_idx)) ...
            - min(Timeline.rawDAQData(:,photodiode_idx)))/2 + ...
            min(Timeline.rawDAQData(:,photodiode_idx));
        % median filter because of weird effect where
        % photodiode dims instead of off for one sample
        % while backlight is turning off
        photodiode_trace = medfilt1(Timeline.rawDAQData(stimScreen_on, ...
            photodiode_idx),10) > photodiode_thresh;
        photodiode_flip = find((~photodiode_trace(1:end-1) & photodiode_trace(2:end)) | ...
            (photodiode_trace(1:end-1) & ~photodiode_trace(2:end)))+1;
        
        photodiode = struct('timestamps',[],'values',[]);
        photodiode.timestamps = stimScreen_on_t(photodiode_flip)';
        photodiode.values = photodiode_trace(photodiode_flip);
end

photodiode_offsets = photodiode.timestamps(photodiode.values == 0);
photodiode_onsets = photodiode.timestamps(photodiode.values == 1);

% Get specific stim onsets by time between last offset and new onset
% (occasionally there a bad frame so flip but not new stim)
refresh_rate_cutoff = 1/5;
stimOn_times = photodiode_onsets( ...
    [1;find(photodiode_onsets(2:end) - photodiode_offsets(1:end-1) > refresh_rate_cutoff) + 1]);

if length(stimOn_times) ~= numel(Protocol.seqnums)
    error('MPEP/Photodiode error: photodiode doesn''t match stim')
end

stimIDs = zeros(size(stimOn_times));
for q = 1:size(Protocol.seqnums,1)
    stimIDs(Protocol.seqnums(q,:)) = q;
end


%% sanity check figures
% stackedplot(Timeline.rawDAQTimestamps, Timeline.rawDAQData(:,[2 3 7 11 13]),...
%     'displaylabels',{'expose clock','pco exposure','blueLED','stimscreen','photodiode'});

ax(1) = subplot(211);
plot(Timeline.rawDAQTimestamps, Timeline.rawDAQData(:,photodiode_idx), ...
    stimScreen_on_t, Timeline.rawDAQData(stimScreen_on, photodiode_idx),...
    stimScreen_on_t, photodiodeOn_filt)
line(Timeline.rawDAQTimestamps([1 end]), [photodiode_thresh photodiode_thresh],'color','k')
legend('raw PD', 'raw PD resample','filtered PD resample', 'th');


ax(2) = subplot(212);
plot(stimScreen_on_t, photodiode_trace, '.-',...
    photodiode.timestamps, photodiode.values, '*'); %....looks like this is miraculuously correct
legend('PD trace', 'PD flip');
ylim([-0.2 1.2])
linkaxes(ax,'x');


%% Load imaging data
if analyseImage
    if cam_color_n == 1
        
        data_path = ['\\zubjects.cortexlab.net\Subjects\' mpep_animal filesep day];
        experiment_path = [data_path filesep num2str(experiment)];
        
        frame_t = readNPY([experiment_path filesep 'svdTemporalComponents_cam2.timestamps.npy']);
        U = readUfromNPY([data_path filesep 'svdSpatialComponents_cam2.npy']);
        V = readVfromNPY([experiment_path filesep 'svdTemporalComponents_cam2.npy']);
        
        framerate = 1./nanmedian(diff(frame_t));
        
        % Detrend data, remove heartbeat if that's within frequency range
        if framerate > 28
            fV = detrendAndFilt(V, framerate);
        else
            highpassCutoff = 0.01; % Hz
            [b100s, a100s] = butter(2, highpassCutoff/(framerate/2), 'high');
            
            dV = detrend(V', 'linear')';
            fV = filter(b100s,a100s,dV,[],2);
        end
        
        avg_im = readNPY([data_path filesep 'meanImage_cam2.npy']);
        
    elseif cam_color_n == 2
        
        % Load in all things as neural (n) or hemodynamic (h)
        
        data_path = ['\\zubjects.cortexlab.net\Subjects\' mpep_animal filesep day];
        experiment_path = [data_path filesep experiment];
        
        tn = readNPY([experiment_path filesep 'svdTemporalComponents_' cam_color_signal '.timestamps.npy']);
        Un = readUfromNPY([data_path filesep 'svdSpatialComponents_' cam_color_signal '.npy']);
        Vn = readVfromNPY([experiment_path filesep 'svdTemporalComponents_' cam_color_signal '.npy']);
        avg_im_n = readNPY([data_path filesep 'meanImage_' cam_color_signal '.npy']);
        
        th = readNPY([experiment_path filesep 'svdTemporalComponents_' cam_color_hemo '.timestamps.npy']);
        Uh = readUfromNPY([data_path filesep 'svdSpatialComponents_' cam_color_hemo '.npy']);
        Vh = readVfromNPY([experiment_path filesep 'svdTemporalComponents_' cam_color_hemo '.npy']);
        avg_im_h = readNPY([data_path filesep 'meanImage_' cam_color_hemo '.npy']);
        
        nSV = size(Vn,1);
        framerate = 1./nanmedian(diff(tn));
        if framerate > 28
            fVn = detrendAndFilt(Vn, framerate);
            fVh = detrendAndFilt(Vh, framerate);
        else
            highpassCutoff = 0.01; % Hz
            [b100s, a100s] = butter(2, highpassCutoff/(framerate/2), 'high');
            
            dVn = detrend(Vn', 'linear')';
            fVn = filter(b100s,a100s,dVn,[],2);
            
            dVh = detrend(Vh', 'linear')';
            fVh = filter(b100s,a100s,dVh,[],2);
        end
        
        % Correct hemodynamic signal in blue from green
        % First need to shift alternating signals to be temporally aligned
        % (shifts neural to hemo)
        % Eliminate odd frames out
        min_frames = min(size(fVn,2),size(fVh,2));
        fVn = fVn(:,1:min_frames);
        fVh = fVh(:,1:min_frames);
        
        fVs_th = SubSampleShift(fVn,1,2);
        
        fVh_Un = ChangeU(Uh,fVh,Un);
        
        hemo_freq = [0.2,3]; %[7 12];
        fVn_hemo = HemoCorrectLocal(Un,fVs_th,fVh_Un,framerate,hemo_freq,3);
        
        % set final U/V to use
        fV = fVn_hemo;
        U = Un;
        avg_im = avg_im_n;
        frame_t = th; % shifted to use hemo color times
        
        attenuation = 10;
        fmax = 7;
        fV = lpFilt_DS(fV, framerate,fmax,attenuation);%23/11/18
        %movieWithTracesSVD(U, fV, frame_t);
    end
    
    
    %% Sparse noise retinotopy (bootstrap across trials to get less noise)
    [Uy,Ux,nSV] = size(U);
    myScreenInfo.windowPtr = NaN; % so we can call the stimulus generation and it won't try to display anything
    stimNum = 1;
    ss = eval([Protocol.xfile(1:end-2) '(myScreenInfo, Protocol.pars(:,stimNum));']);
    stim_screen = cat(3,ss.ImageTextures{:});
    ny = size(stim_screen,1);
    nx = size(stim_screen,2);
    switch lower(photodiode_type)
        case 'flicker'
            % Check for case of mismatch between photodiode and stimuli:
            % odd number of stimuli, but one extra photodiode flip to come back down
            if mod(size(stim_screen,3),2) == 1 && ...
                    length(photodiode.timestamps) == size(stim_screen,3) + 1
                photodiode.timestamps(end) = [];
                photodiode.values(end) = [];
                warning('Odd number of stimuli, removed last photodiode');
            end
            
            % If there's still a mismatch, break
            if size(stim_screen,3) ~= length(photodiode.timestamps)
                warning([num2str(size(stim_screen,3)) ' stimuli, ', ...
                    num2str(length(photodiode.timestamps)) ' photodiode pulses']);
                
                % Try to estimate which stim were missed by time difference
                photodiode_diff = diff(photodiode.timestamps);
                max_regular_diff_time = prctile(diff(photodiode.timestamps),99);
                skip_cutoff = max_regular_diff_time*2;
                photodiode_skip = find(photodiode_diff > skip_cutoff);
                est_n_pulse_skip = ceil(photodiode_diff(photodiode_skip)/max_regular_diff_time)-1;
                stim_skip = cell2mat(arrayfun(@(x) photodiode_skip(x):photodiode_skip(x)+est_n_pulse_skip(x)-1, ...
                    1:length(photodiode_skip),'uni',false));
                
                if isempty(est_n_pulse_skip) || length(photodiode.timestamps) + sum(est_n_pulse_skip) ~= size(stim_screen,3)
                    error('Can''t match photodiode events to stimuli')
                end
            end
            
            stim_times = photodiode.timestamps;
            
        case 'steady'
            % If the photodiode is on steady: extrapolate the stim times
            if length(photodiode.timestamps) ~= 2
                error('Steady photodiode, but not 2 flips')
            end
            stim_duration = diff(photodiode.timestamps)/size(stim_screen,3);
            stim_times = linspace(photodiode.timestamps(1), ...
                photodiode.timestamps(2)-stim_duration,size(stim_screen,3))';
            
    end
    
    % Get average response to each stimulus
    surround_window = [0,0.5]; % 6s = [0.1,0.5], 6f = [0.05,0.2]
    framerate = 1./nanmedian(diff(frame_t));
    surround_samplerate = 1/(framerate*1);
    surround_time = surround_window(1):surround_samplerate:surround_window(2);
    response_n = nan(ny,nx);
    response_grid = cell(ny,nx); %V
    response_grid_r = cell(ny,nx); %real space
    for px_y = 1:ny
        for px_x = 1:nx
            
            % Use first frame of dark or light stim
            align_stims = (stim_screen(px_y,px_x,2:end)~= 0) & ...
                (diff(stim_screen(px_y,px_x,:),[],3) ~= 0);
            align_times = stim_times(find(align_stims)+1);
            %         keyboard
            align_times = align_times(round(length(align_times)/2):end);
            
            response_n(px_y,px_x) = length(align_times);
            %         keyboard
            % Don't use times that fall outside of imaging
            align_times(align_times + surround_time(1) < frame_t(2) | ...
                align_times + surround_time(2) > frame_t(end)) = [];
            %         keyboard
            % Get stim-aligned responses, 2 choices:
            
            % 1) Interpolate times (slow - but supersamples so better)
            %         align_surround_times = bsxfun(@plus, align_times, surround_time);
            %         peri_stim_v = permute(mean(interp1(frame_t,fV',align_surround_times),1),[3,2,1]);
            
            % 2) Use closest frames to times (much faster - not different)
            align_surround_times = bsxfun(@plus, align_times, surround_time);
            frame_edges = [frame_t,frame_t(end)+1/framerate];
            align_frames = discretize(align_surround_times,frame_edges);
            
            align_frames(any(isnan(align_frames),2),:) = [];
            %         keyboard
            % Define the peri-stim V's as subtracting first frame (baseline)
            peri_stim_v = bsxfun(@minus, ...
                reshape(fV(:,align_frames)',size(align_frames,1),size(align_frames,2),[]), ...
                reshape(fV(:,align_frames(:,1))',size(align_frames(:,1),1),size(align_frames(:,1),2),[]));
            % trial x time x nSVD
            
            %         keyboard
            mean_peri_stim_v = permute(mean(peri_stim_v,2),[3,1,2]); %mean across time. [nSVD x trial]
            %         keyboard
            % Save V's
            response_grid{px_y,px_x} = mean_peri_stim_v;
            
            %                 eventLabels = ones(1, length(align_times));
            %  pixelTuningCurveViewerSVD(Ud, fV, frame_t, align_times, eventLabels, surround_window)
            
            %         %% response in real space
            %         peri_stim_v_c = permute(peri_stim_v, [3 1 2]);
            %         peri_stim_r_c = svdFrameReconstruct(U, reshape(peri_stim_v_c, size(peri_stim_v_c, 1), []));
            %         dperi_stim_r = max(diff(reshape(peri_stim_r_c, size(peri_stim_r_c,1),size(peri_stim_r_c,2), size(peri_stim_v_c,2), size(peri_stim_v_c,3)),3),0); %diff & rectify
            %
            %         response_grid_r{px_y,px_x} = mean(dperi_stim_r, 4);
            %
            %             subplot(ny,nx,px_x+nx*(px_y-1));
            %             imagesc(mean(response_grid_r{px_y,px_x},3));
            %             title(['x:' num2str(px_x) ', y:' num2str(px_y)])
            %             axis image off;
            
        end
    end
    
    
    
    % Get position preference for every pixel
    U_downsample_factor = 1; %2 if max method
    screen_resize_scale = 1; %3 if max method
    filter_sigma = (screen_resize_scale*2);
    
    % Downsample U
    use_u_y = 1:Uy;
    Ud = imresize(U(use_u_y,:,:),1/U_downsample_factor,'bilinear');
    
    % Convert V responses to pixel responses
    use_svs = 1:size(U,3); %number of svd components
    n_boot = 10;
    
    response_mean_boostrap = cellfun(@(x) bootstrp(n_boot,@mean,x')',response_grid,'uni',false);
    %response_mean_boostrap_r = cellfun(@(x) bootstrp(n_boot,@mean,x')',response_grid_r,'uni',false);
    
    % (to split trials instead of bootstrap)
    %split_trials = cellfun(@(x) shake(discretize(1:size(x,2),round(linspace(1,size(x,2),n_boot+1)))),response_grid,'uni',false);
    %response_mean_boostrap = cellfun(@(x,y) grpstats(x',y','mean')',response_grid,split_trials,'uni',false);
    use_method = 'max'; % max or com
    vfs_boot = nan(size(Ud,1),size(Ud,2),n_boot);
    for curr_boot = 1:n_boot
        %     keyboard
        response_mean = cell2mat(cellfun(@(x) x(:,curr_boot),response_mean_boostrap(:),'uni',false)');
        stim_im_px = reshape(permute(svdFrameReconstruct(Ud(:,:,use_svs),response_mean(use_svs,:)),[3,1,2]),ny,nx,[]);
        gauss_filt = fspecial('gaussian',[ny,nx],filter_sigma);
        stim_im_smoothed = imfilter(imresize(stim_im_px,screen_resize_scale,'bilinear'),gauss_filt);
        
        switch use_method
            case 'max'
                % Upsample each pixel's response map and find maximum
                [~,mi] = max(reshape(stim_im_smoothed,[],size(stim_im_px,3)),[],1);
                [m_y,m_x] = ind2sub(size(stim_im_smoothed),mi);
                m_yr = reshape(m_y,size(Ud,1),size(Ud,2));
                m_xr = reshape(m_x,size(Ud,1),size(Ud,2));
                
            case 'com'
                % Conversely, do COM on original^2
                [xx,yy] = meshgrid(1:size(stim_im_smoothed,2),1:size(stim_im_smoothed,1));
                m_xr = reshape(sum(sum(bsxfun(@times,stim_im_smoothed.^2,xx),1),2)./sum(sum(stim_im_smoothed.^2,1),2),size(Ud,1),size(Ud,2));
                m_yr = reshape(sum(sum(bsxfun(@times,stim_im_smoothed.^2,yy),1),2)./sum(sum(stim_im_smoothed.^2,1),2),size(Ud,1),size(Ud,2));
        end
        
        % Calculate and plot sign map (dot product between horz & vert gradient)
        
        % 1) get gradient direction
        [~,Vdir] = imgradient(imgaussfilt(m_yr,1));
        [~,Hdir] = imgradient(imgaussfilt(m_xr,1));
        
        % 3) get sin(difference in direction) if retinotopic, H/V should be
        % orthogonal, so the closer the orthogonal the better (and get sign)
        angle_diff = sind(Vdir-Hdir);
        angle_diff(isnan(angle_diff)) = 0;
        
        vfs_boot(:,:,curr_boot) = angle_diff;
        
        disp(curr_boot);
        
    end
    
    vfs_median = imgaussfilt(nanmedian(vfs_boot,3),2);
    
    figure('Name',mpep_animal);
    ax1 = axes;
    subplot(1,2,1,ax1);
    imagesc(vfs_median);
    caxis([-1,1]);
    axes(ax1); axis image off;
    colormap(colormap_BlueWhiteRed)
    mcolorbar(gca, .5);
    
    ax2 = axes;
    subplot(1,2,2,ax2);
    h1 = imagesc(ax2,avg_im_h(use_u_y,:)); %avg_im_n
    caxis(ax2,[0 prctile(avg_im_h(:),98)]);hold on
    set(ax2,'Visible','off');
    axes(ax2); axis image off;
    
    ax3 = axes;
    subplot(1,2,2,ax3);
    h2 = imagesc(ax3,vfs_median);
    colormap(ax3,colormap_BlueWhiteRed);
    caxis([-1,1]);
    set(ax3,'Visible','off');
    axes(ax3); axis image off;
    set(h2,'AlphaData',mat2gray(abs(vfs_median))*0.6); %0.5
    
    colormap(ax2,gray);
    
    savePaperFigure(gcf, ['vfs ' mpep_animal ' ' day ' ' experiment]);
end

