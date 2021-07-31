function expt = grabStimTimes2p(expt, getAllPhdFlips, TSDir, ISI, bklightCtrl, thr, detectTdiff)
%times = getStimTimesWF
% returns times of visual stimulus onset from Thorsync data
% 
% This function is backward compatible: if "syncSquare" and "stimScreen"
% are not recorded, use only "Screen"
%
% INPUT:
% expt.subject, expt.expDate, expt.expNum
% getAllPhdFlips: if 1, return times.frameTimes
% ISI: interval between end of previous trial to the beginning of the next trial. used to detect stimulus onset
% NOTE: this definition is different to that used in grabStimTimes2ps
% bklightCtrl: whether the backlight is controled by arduino or not
% (default: true).
%
% OUTPUT:
% expt.stimTimes.onset
% expt.stimTimes.offset
% expt.stimTimes.frameTimes{iStim}
%
% assumes flicker SyncSquare & following thorsync entries:
% Screen ... measured photodiode signal from screen(corresponding to "photoDiode" in wf).
% syncSquare ... copy of syncsquare (without backlight ctrl) (from vs pc)
% stimScreen ... signal sent to backlight control (from arduino)
%
% 10/7/20 DS created from grabStimTimesWF
%

%TODO check if this is compatible when syncsquare = steady

verbose = true;

if nargin < 7
    detectTdiff = 1; 
end
if nargin < 6 || isempty(thr)
    thr = 0.5; %used only for grabStimTimes2ps
end

if nargin < 5 || isempty(bklightCtrl)
    bklightCtrl = 1;
end

if nargin < 4 || isempty(ISI)
    ISI = 1/5; %[s] inter-repeat interval
end

if nargin < 3 || isempty(TSDir)
    %thisDate = expt.expDate(1:10);
    %thisSeries = str2num(expt.expDate(12:end));
    
    TSDir = 'D:\thorimagedata\';
end
if nargin < 2 || isempty(getAllPhdFlips)
    getAllPhdFlips = 0;
end


disp('Grabbing stim times ');
if isfield(expt, 'stimTimes')
    expt = rmfield(expt, 'stimTimes');%15/12/20 for grabStimTimes2ps
end

%% Load thorsync data
fnam = fullfile(TSDir, expt.subject, expt.expDate, num2str(expt.expNum), 'Episode001.h5');
if verbose; disp('Loading thorsync...'); end
lineNames = {'stimScreen','stimScreeCraig','syncSquare','Screen'};
[syncDataOut] = LoadSyncEpisodeFunction(fnam, lineNames);

%%for backward compatibility
if ~isfield(syncDataOut,'stimScreen')
     expt = grabStimTimes2ps( expt, 0, TSDir, thr, ISI );
    return;
end
    
%omit these?:
expt.TSfile = fnam;
expt.TSDir = TSDir;
    
tstime = syncDataOut.time;

%% Get wheter syncSquare on screen was flickering
if isfield(syncDataOut,'syncSquare')
    syncSquare = double(syncDataOut.syncSquare);
    %syncSquare_flicker = max(syncSquare) - min(syncSquare) > 2;
    syncSquare_thresh = max(syncSquare)/2;
    syncSquare_white_on = syncSquare > syncSquare_thresh; %[0 1]
    syncSquare_white_on_t = tstime(syncSquare_white_on);
    
    white_off = find(diff(syncSquare_white_on_t) > 2/55);% assume screen resolution at 60Hz
    syncSquare_white_start_t = [syncSquare_white_on_t(1); syncSquare_white_on_t(white_off + 1)]; % ~stimulus start
    syncSquare_white_end_t = [syncSquare_white_on_t(white_off); syncSquare_white_on_t(end)]; % ~stimulus end
    syncSquare_withinStim = logical(zeros(numel(tstime),1));
    for irepeat = 1:numel(syncSquare_white_start_t)
        syncSquare_withinStim((tstime>=syncSquare_white_start_t(irepeat)) & ...
            (tstime<=syncSquare_white_end_t(irepeat))) = 1;
    end
    syncSquare_withinStim = logical(syncSquare_withinStim);
    syncSquare_black_on = logical(syncSquare_withinStim .* ~syncSquare_white_on);
    syncSquare_black_on_t = tstime(syncSquare_black_on);
end



%% Get whether backlight was flickering
% this is equal to input to backlight control
% this channel is not very useful??
if isfield(syncDataOut, 'stimScreen') && bklightCtrl
    stimScreen = double(syncDataOut.stimScreen);
    %stimScreen_flicker = max(stimScreen) - min(stimScreen) > 2;
    stimScreen_thresh = 0.5;%for 2p
    stimScreen_on = stimScreen > stimScreen_thresh; %[0 1]
else
    stimScreen_on = ones(length(tstime),1);
end



%% Get flicker or steady photodiode
%     photodiode_type = myScreenInfo.SyncSquare.Type;
% THE CODE BELOW SHOULD WORK IN BOTH STEADY AND FLICKER BUT NOT YET CHECKED
phd = syncDataOut.Screen;

if bklightCtrl
    %     switch lower(photodiode_type)
    %         case 'flicker' %assumes black background (RigInfo.BackgroundColor=[0 0 0])
    %when backlight is off. complete darkness
    bkOff = mode(phd(~stimScreen_on));
    
    %when backlight is on and stimulus is NOT presented ...
    bkOn_nostim = mode(phd(logical(~syncSquare_withinStim.*stimScreen_on)));
    
    screenOnTh = ( bkOff + bkOn_nostim)/2;
else
    screenOnTh = 0; %7/7/20
end

%% adjust time difference between syncsquare(sent from vs) and photodiode (measured on screen)
photodiode_on = phd > screenOnTh;
%more robust?: photodiode_on = photodiode_on .* stimScreen_on;

photodiode_thresh = max(phd)/2; %make this lower?
photodiode_white_on = phd > photodiode_thresh; %[0 1]
photodiode_start_t = min(tstime(photodiode_white_on));

%time difference between syncsquare and photodiode
if detectTdiff
    tdiff = min(syncSquare_white_start_t) - photodiode_start_t ;
else
    tdiff = 0;
end
disp(['time difference between measured photodiode and signal from vs is ' num2str(tdiff) '[s]']);
%     syncSquare_white_on_t_adj = syncSquare_white_on_t + tdiff;
%     syncSquare_white_on_adj = tstime
tdiffIdx = round(tdiff/median(diff(tstime)));
if tdiff > 0
    syncSquare_white_on_adj = [syncSquare_white_on(tdiffIdx:end); zeros(tdiffIdx-1,1)];
elseif tdiff < 0
    syncSquare_white_on_adj = [zeros(-tdiffIdx,1); syncSquare_white_on(1:end+tdiffIdx)];
else 
    syncSquare_white_on_adj = syncSquare_white_on; %4/11/20
end
%     %sanity check
%     plot(tstime, syncSquare_white_on)
%     hold on
%     plot(tstime, syncSquare_white_on_adj)
%     plot(tstime, phd)

%time when photodiode detect backlight on
photodiode_on_t = tstime(photodiode_on)';
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
photodiode.timestamps = photodiode_on_t(photodiode_flip)';
photodiode.values = photodiode_trace(photodiode_flip);
%< this is used to detect:
% 1, stimulus on/offset of each repeat
% 2, onset of each stimulus frame, specified in vs


%% sanity check figures
if verbose
    figure;
    ax(1) = subplot(311);
    plot(tstime, phd, photodiode_on_t, phd(photodiode_on));
    grid on;
    line(tstime([1 end]), [screenOnTh screenOnTh],'color','k')
    legend('raw PD', 'raw PD resample','th');
    
    ax(2) = subplot(312);
    plot(tstime,  syncSquare, tstime,  syncSquare_white_on_adj);
    grid on;
    legend('syncSquare (before adj)','syncSquare (after adj)');
    
    ax(3) = subplot(313);
    plot(photodiode_on_t, photodiode_trace, '.-',...
        photodiode.timestamps, photodiode.values, '*'); %....looks like this is miraculuously correct
    grid on;
    legend('PD trace', 'PD flip');
    ylim([-0.2 1.2])
    
    linkaxes(ax,'x');
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
stimOn_idx = [1;find(photodiode_onsets(2:end) - photodiode_offsets(1:end-1) > ISI) + 1];
stimOff_idx = [stimOn_idx(2:end)-1; length(photodiode_offsets)];
stimOn_times = photodiode_onsets(stimOn_idx);
stimOff_times = photodiode_offsets(stimOff_idx);
if verbose
    ax(3) = subplot(313);
    hold on;
    plot(stimOn_times, 0, 'o', stimOff_times, 0, 'x');
    linkaxes(ax,'x');
end


times.onset = stimOn_times;
times.offset = stimOff_times; %maynot be accurate
if getAllPhdFlips
    for iStim = 1:length(stimOn_times)
        tmpUp = photodiode_onsets(photodiode_onsets>=times.onset(iStim) ...
            & photodiode_onsets<=times.offset(iStim));
        tmpDown = photodiode_offsets(photodiode_offsets>=times.onset(iStim) ...
            & photodiode_offsets<=times.offset(iStim));
        times.frameTimes{iStim,1} = sort([tmpUp(:); tmpDown(:)], 'ascend');
        
    end
end

expt.stimTimes = times;