function expt = grabStimTimesWF(expt, getAllPhdFlips, TLDir, ISI, bklightCtrl, verbose)
%times = grabStimTimesWF(expt, getAllPhdFlips, TLDir, ISI, bklightCtrl, verbose)
%returns times of visual stimulus onset from Timeline data.
%note this function does not read protocol file, so the order of output
%variables are not in the order stored in p-file
%
% INPUT:
% expt.subject, expt.expDate, expt.expNum
% getAllPhdFlips: if 0,don't return times.frameTimes, otherwise
% getAllPhdFlips(1) = nstims, getAllPhdFlips(2) = nrepeats, used for times.frameTimes
% ISI: inter-stimulus interval[s]. used to detect stimulus onset
% bklightCtrl: whether the backlight is controled by arduino or not
% (default: true).
%
% OUTPUT:
% expt.stimTimes.onset: start time of each trial
% expt.stimTimes.offset: stop time of each trial
% expt.stimTimes.frameTimes{iStim}
%
% assumes flicker SyncSquare & following timeline entries:
% photoDiode ... measured photodiode signal from screen.
% syncSquare ... copy of syncsquare (without backlight ctrl) (from vs pc)
% stimScreen ... signal sent to backlight control (from arduino)
%
% 11/6/20 DS created from lilrig_load_experiment.m and getStimTimesWF
% 7/7/20 added bklightCtrl
% 7/7/20 fixed time difference between syncSquare and photodiode
% 10/7/20 definition of TLDir changed
% 15/10/20 fixed bug when #photodiode onset = #photodiode offset + 1
% 10/11/20 detect tdiff in every trial


%% TODO:
% detect Tdiff for each trial
% judge whether backlight is controlled with frequency analogy of
% photodiode and camera exposure (to know camera sampling freq)

detectTDiff = 1; %detect time difference between photodiode and syncsquare 

if nargin < 6 || isempty(verbose)
    verbose = 0; % Display progress or not
end

if nargin < 5 || isempty(bklightCtrl)
    bklightCtrl = 1;
end

if nargin < 4 || isempty(ISI)
    ISI = 1/5; %[s] inter-repeat interval
end

if nargin < 3 || isempty(TLDir)
    %thisDate = expt.expDate(1:10);
    %thisSeries = str2num(expt.expDate(12:end));
    %TLDir = fileparts(dat.expFilePath(expt.subject, thisDate, thisSeries, expt.expNum, 'timeline','master'));
    TLDir = '\\storage.erc.monash.edu.au\shares\MNHS-dshi0006\Subjects';
end
if nargin < 2
    getAllPhdFlips = 0;
end




%% Load timeline

% [timeline_filename, timeline_exists] = AP_cortexlab_filename(expt.subject,expt.expDate,expt.expNum,'timeline');
% TLDir = fullfile('C:\Users\Experiment\Documents\MATLAB\Data\Subjects\',...
%     expt.subject,expt.expDate,num2str(expt.expNum));
filename_TL = sprintf('%s_%d_%s_Timeline.mat',expt.expDate, expt.expNum, expt.subject);
timeline_filename = fullfile(TLDir, expt.subject, expt.expDate, num2str(expt.expNum), filename_TL);

if ~exist(timeline_filename,'file')
    error([expt.subject ' ' expt.expDate ': no timeline']);
else
    if verbose; disp('Loading timeline...'); end
    
    load(timeline_filename);
    
    tltime = Timeline.rawDAQTimestamps';
    
    %% Get wheter syncSquare on screen was flickering
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
    if any(stimScreen_idx) && bklightCtrl
        stimScreen_flicker = max(Timeline.rawDAQData(:,stimScreen_idx)) - ...
            min(Timeline.rawDAQData(:,stimScreen_idx)) > 2;
        %stimScreen_thresh = max(Timeline.rawDAQData(:,stimScreen_idx))/2;
        stimScreen_thresh = 2.5;%this way compatible with backlight is always on
        stimScreen_on = Timeline.rawDAQData(:,stimScreen_idx) > stimScreen_thresh; %[0 1]
    else
        stimScreen_on = ones(Timeline.rawDAQSampleCount,1);
    end
    
    
    
    %% Get flicker or steady photodiode
    
    
    % detect camera sampling rate = backlight freq
    %     L = length(ROItrace);
    %     NFFT = 2^nextpow2(L);
    %     [Pxx,F] = pwelch(ROItrace-mode(ROItrace),[],[],NFFT,Fs);
    %     %    plot(F, 10*log10(Pxx));
    %     theseFreqInds = find(F>1);
    %     [~, maxInds] = max(Pxx(theseFreqInds));
    %     hemoFreq = [max(1,F(theseFreqInds(maxInds))-2) min(F(theseFreqInds(maxInds))+2, 0.99*Fs/2)];%15/10/20

    % judge if there is frequency component in photodiode signal
    
    
    %     photodiode_type = myScreenInfo.SyncSquare.Type;
    % THE CODE BELOW SHOULD WORK IN BOTH STEADY AND FLICKER BUT NOT YET CHECKED
    photodiode_idx = strcmp({Timeline.hw.inputs.name}, 'photoDiode');
     
    if bklightCtrl
        %     switch lower(photodiode_type)
        %         case 'flicker' %assumes black background (RigInfo.BackgroundColor=[0 0 0])
        %when backlight is off. complete darkness
        bkOff = mode(Timeline.rawDAQData(~stimScreen_on,photodiode_idx));
        
        %when backlight is on and stimulus is NOT presented ...
        bkOn_nostim = mode(Timeline.rawDAQData(logical(~syncSquare_withinStim.*stimScreen_on),photodiode_idx));
        
        screenOnTh = ( bkOff + bkOn_nostim)/2;
    else
        screenOnTh = 0; %25/9/20 this works better than one line below
        %screenOnTh = mean(Timeline.rawDAQData(:,photodiode_idx)); %25/9/20
    end
    
    %% adjust time difference between syncsquare(sent from vs) and photodiode (measured on screen)
    % 10/11/20 detect tdiff in every trial
    phd = Timeline.rawDAQData(:,photodiode_idx)';
    photodiode_on = phd > screenOnTh;
    %more robust?: photodiode_on = photodiode_on .* stimScreen_on;
    srate = 1/median(diff(tltime));
    
    photodiode_thresh = max(phd)/2; %make this lower?
    photodiode_white_on = phd > photodiode_thresh; %[0 1]
    photodiode_white_on_t = tltime(photodiode_white_on);
    
    syncSquare_white_on_adj = zeros(length(tltime),1);
    for iStart = 1:length(syncSquare_white_start_t)
        syncStart_t = syncSquare_white_start_t(iStart);
        margin_t = [syncStart_t - ISI/2 syncStart_t + ISI/2];
        %margin_idx = intersect(find(tltime > margin_t(1)), find(tltime < margin_t(2)));
        pdStart_idx = min(intersect(find( photodiode_white_on_t > margin_t(1)), find(photodiode_white_on_t < margin_t(2))));
        pdStart_t = photodiode_white_on_t(pdStart_idx);
        
        tidx = find(tltime == syncStart_t) : find(tltime == syncSquare_white_end_t(iStart));
         
        tidx_margin = tidx(1) - round(0.5*ISI*srate): tidx(end) + round(0.5*ISI*srate);
        
        %time difference between syncsquare and photodiode
        if detectTDiff
            tdiff = syncStart_t - pdStart_t; %this is not robust
            
            %             %TODO: compute xcorr to estimate delay ... so far its not
            %             %working better than detecting onsets
            %             [c, lags] = xcorr(syncSquare_white_on(tidx_margin)-mean(syncSquare_white_on(tidx_margin)), ...
            %                 photodiode_white_on(tidx_margin)-mean(photodiode_white_on(tidx_margin)));
            %             [~,lagidx] = max(c);
            %             tdiff = lags(lagidx)/srate;
            
        else
            tdiff = 0; %17/9/20
        end
        disp(['Repeat:' num2str(iStart) ', time difference between measured photodiode and signal from vs is ' num2str(tdiff) '[s]']);
        %     syncSquare_white_on_t_adj = syncSquare_white_on_t + tdiff;
        %     syncSquare_white_on_adj = tltime
        tdiffIdx = round(tdiff*srate);
        tidx_adjusted = tidx - tdiffIdx;
        syncSquare_white_on_adj(tidx_adjusted) =  syncSquare_white_on(tidx);
        
        %         if tdiff > 0
        %             syncSquare_white_on_adj = [syncSquare_white_on(tdiffIdx:end); zeros(tdiffIdx-1,1)];
        %         elseif tdiff < 0
        %             syncSquare_white_on_adj = [zeros(-tdiffIdx,1); syncSquare_white_on(1:end+tdiffIdx)];
        %         elseif tdiff == 0
        %             syncSquare_white_on_adj = syncSquare_white_on;
        %         end
        %     %sanity check
        %     plot(tltime, syncSquare_white_on)
        %     hold on
        %     plot(tltime, syncSquare_white_on_adj)
        %     plot(tltime, phd)
    end
    
    %time when photodiode detect backlight on
    photodiode_on_t = tltime(photodiode_on);
    %             photodiodeOn_filt = Timeline.rawDAQData(photodiode_on, photodiode_idx); %medfilt1(Timeline.rawDAQData(stimScreen_on, photodiode_idx),5);
    %             photodiode_thresh = median(photodiodeOn_filt);
    
    % 1 when flicker syncsquare is HIGH, 0 when flicker syncsquare is LOW
    %             photodiode_trace = photodiodeOn_filt > photodiode_thresh;
    % photodiode_trace = syncSquare_white_on(photodiode_on); %4/5 use syncSquare w/b!
    photodiode_trace = syncSquare_white_on_adj(photodiode_on); %7/7/20
    
    %this judges whether the photodiode signal flipped or continued
    photodiode_flip = find((~photodiode_trace(1:end-1) & photodiode_trace(2:end)) | ...
        (photodiode_trace(1:end-1) & ~photodiode_trace(2:end)))+1; %time idx when signal flips 0 to 1 or 1 to 0 ??
    
    
    photodiode = struct('timestamps',[],'values',[]);
    photodiode.timestamps = photodiode_on_t(photodiode_flip);
    photodiode.values = photodiode_trace(photodiode_flip)';
    %< this is used to detect:
    % 1, stimulus on/offset of each repeat
    % 2, onset of each stimulus frame, specified in vs
    
    
    %% sanity check figures
    if verbose
        figure;
        ax(1) = subplot(311);
        plot(tltime, phd, photodiode_on_t, phd(photodiode_on));
        grid on;
        line(tltime([1 end]), [screenOnTh screenOnTh],'color','k')
        legend('raw PD', 'raw PD resample','th');
        
        ax(2) = subplot(312);
        plot(tltime,  Timeline.rawDAQData(:,syncSquare_idx),...
            tltime,  syncSquare_white_on_adj);
        grid on;
        legend('syncSquare (before adj)','syncSquare (after adj)');
        
        ax(3) = subplot(313);
        plot(photodiode_on_t, photodiode_trace, '.-',...
            photodiode.timestamps, photodiode.values, '*'); %....looks like this is miraculuously correct
        grid on;
        legend('PD trace', 'PD flip');
        ylim([-0.2 1.2]);
        linkaxes(ax,'x');
    end
    
    %         case 'steady'
    %             warning('dont use Steady syncSquare!!');
    %     end
    photodiode_offsets = photodiode.timestamps(photodiode.values == 0); %maybe missing/adding last frame of each repeat
    photodiode_onsets = photodiode.timestamps(photodiode.values == 1);   
    
    %heuristics numbers of onsets and offsets do not match 14/10/20
    if (numel(photodiode_onsets) == numel(photodiode_offsets) + 1)
    upDuration = photodiode_offsets(1:end)-photodiode_onsets(1:end-1);
        if isempty(find(upDuration < 0, 1))
            photodiode_offsets = [photodiode_offsets; photodiode_onsets(end)+median(upDuration)];
        end
    end
    
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
    %     stimOn_idx = [1; find(photodiode_onsets(2:end) - photodiode_offsets(1:end-1) > ISI) + 1];
    stimOn_idx = [1; find(diff(photodiode_onsets)>ISI) + 1]; %15/1/20 more robust than the above if long (200s) trial
    
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
if getAllPhdFlips ~= 0
    nstims = getAllPhdFlips(1);
    nrepeats = getAllPhdFlips(2);
    for iRepeat = 1:nrepeats
        for iStim = 1:nstims%length(stimOn_times)
            
            idx = nstims * (iRepeat - 1) + iStim; %9/11/20
            
            tmpUp = photodiode_onsets(photodiode_onsets>=times.onset(idx) ...
                & photodiode_onsets<=times.offset(idx));
            tmpDown = photodiode_offsets(photodiode_offsets>=times.onset(idx) ...
                & photodiode_offsets<=times.offset(idx));
            times.frameTimes{iStim, iRepeat} = sort([tmpUp(:); tmpDown(:)], 'ascend')';
        end
    end
end

expt.stimTimes = times;