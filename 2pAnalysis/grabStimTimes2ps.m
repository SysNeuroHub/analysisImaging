function expt = grabStimTimes2ps( expt, getAllPhdFlips, TSDir, thr, ISI )
%expt = grabStimTimes2ps( expt, pdFlag, getAllPhdFlips, TSDir, thr, ISI)

%13/1/20 created from grabStimTimes for Timeline version
%10/7/20 renamed from grabStimTimes
%
%this version only detects onset of each repeat insteady pd
%not compatible with flicker syncSquare
%hence getAllPhdFlips must be false 

subject = expt.subject;
expDate = expt.expDate;
expr = expt.expNum; %can be a vector

if nargin < 5
    ISI = 1;
end

if nargin < 4 
    thr = 0.5;
end
if nargin < 3
    TSDir = 'D:\thorimagedata\'; %where ThorSync data is saved or copied
end

if nargin < 2
    getAllPhdFlips = false;
end

disp('Retrieving stimulus presentation times from ThorSync...')
for iexp = 1:length(expt.expNum)
    
    %timelineFile = strcat(expDate,'_',num2str(expr(iexp)),'_',subject,'_','Timeline.mat');
    fnam = fullfile(TSDir, subject, expDate, num2str(expr), 'Episode001.h5');
    [syncDataOut] = LoadSyncEpisodeFunction(fnam);
    
    expt.TSfile{iexp} = fnam;
    expt.TSDir = TSDir;
    
    
    %ind = strcmp({Timeline.hw.inputs.name}, 'photoDiode');
    time = syncDataOut.time;
    try
        phd = syncDataOut.Screen;
    catch err
        phd = syncDataOut.Screen1; %20/11/20
    end
    %phd = [fliplr(phd); phd; fliplr(phd)];
    
    
    %% resample to eliminate 12kHz
    fs_rs = 1e3;
    [phd_rs, time_rs] = resample(phd, time, fs_rs);
    
    phdMax = prctile(phd_rs,99);%max(phd);
    phdMin = prctile(phd_rs,1);%min(phd);
    
    phd_rs = (phd_rs-phdMin)./(phdMax-phdMin);
    
    
    %% notch filtering of 50Hz
    wo = 50/(fs_rs/2);
    bw = 300*wo/35;
     [b,a] = iirnotch(wo,bw);
%             h=fvtool(b,a);
%             h.Fs = fs_rs;
%             h.FrequencyScale= 'Log';
    phd_f = filtfilt(b,a, phd_rs);
    
    %further resample to eliminate 50Hz
    fs_f = 50;
    [phd_f, time_f] = resample(phd_f, time_rs, fs_f);
    
    phdMax = prctile(phd_f,99);%max(phd);
    phdMin = prctile(phd_f,1);%min(phd);
    
    phd_f = (phd_f-phdMin)./(phdMax-phdMin);
    
    
    %thr = 0.35; % Set this lower as in some cases there are slow changes in phd amplitude. Check this if the number of frames seems incorrect
    %thr = 0.7; %when syncSquare includes grey
    
    if strcmp(thr, 'userInput')
        figure;
        plot(time_rs, phd_rs,time_f, phd_f);grid minor;
        thr = input('What is the threshold for photodiode screen?\n');
        close;
    end
    
    above = phd_rs>thr;
    deltas = [0; diff(above)];
    
    above_f = phd_f>thr;
    deltas_f = [0; diff(above_f)];
    
    upIdx = find(deltas==1); %high time res but includes errorneous events
    upIdx_f = find(deltas_f==1);%low time res but does not include errorneous event
    goingUpTimes = [];
    %goingUpTimes = time_rs(deltas==1);
    dt = 0.25; %search the up time within this time range
    for uu = 1:length(upIdx_f)
        trange = [time_f(upIdx_f(uu))-dt time_f(upIdx_f(uu))+dt];
        inTrange = find((time_rs(upIdx)>trange(1)) .* (time_rs(upIdx)<trange(2)));
        if ~isempty(inTrange)
            [~,bestone] = min(abs(upIdx(inTrange) - upIdx_f(uu)));
            goingUpTimes = [goingUpTimes; time_rs(upIdx(inTrange(bestone)))];
        end
    end
    
    
    
    downIdx = find(deltas == -1); %high time res but includes errorneous events
    downIdx_f = find(deltas_f == -1);%low time res but does not include errorneous event
    goingDownTimes = [];
    %goingDownTimes = time_rs(deltas==-1);
    for uu = 1:length(downIdx_f)
        trange = [time_f(downIdx_f(uu))-dt time_f(downIdx_f(uu))+dt];
        inTrange = find((time_rs(downIdx)>trange(1)) .* (time_rs(downIdx)<trange(2)));
        if ~isempty(inTrange)
            [~,bestone] = min(abs(downIdx(inTrange) - downIdx_f(uu)));
            goingDownTimes = [goingDownTimes; time_rs(downIdx(inTrange(bestone)))];
        end
    end
    
    %% trim too close events
    % while 1
    diffUpTimes = diff(goingUpTimes);
    okIdx = [1; 1 + find(diffUpTimes > ISI)];
    goingUpTimes = goingUpTimes(okIdx);
    %     if isempty(diff(goingUpTimes) <= ISI)
    %         continue;
    %     end
    % end
    
    times.onset = goingUpTimes;
    % for ii = 1:length(diffUpTimes)
    %     if diffUpTimes(ii) > 0.1 %[s]
    %         time.onset = [times.onset; goingUpTimes(ii)];
    %     end
    % end
    
    
    %     times.onset=[];
    %     for iUDP=1:length(startIdx)
    %         tmp=min(goingUpTimes(goingUpTimes>Timeline.mpepUDPTimes(startIdx(iUDP))));
    %         times.onset=[times.onset; tmp];
    %     end
    
    %     % the last DOWN before the StimEnd UDP
    %     times.offset=[];
    %     delayConst = 0.05; % a delay constant, for the case when the 'StimEnd' udp arrived before the stimulus finished to play
    %     for iUDP=1:length(endIdx)
    %         tmp=max(goingDownTimes(goingDownTimes<(Timeline.mpepUDPTimes(endIdx(iUDP)) + delayConst)));
    %         times.offset=[times.offset; tmp];
    %     end
    diffDownTimes = diff(goingDownTimes);
    okIdx = [1; 1 + find(diffDownTimes > ISI)];
    goingDownTimes = goingDownTimes(okIdx);
    
    % offset must be AFTER onset %10/7/20
    times.offset = nan(length(times.onset),1);
    for tt = 1:length(times.onset)
        times.offset(tt) = goingDownTimes(min(find(goingDownTimes - times.onset(tt) > 0)));
    end
    %times.offset = goingDownTimes;
    
    if getAllPhdFlips
        for iStim = 1:length(times.onset)
            tmpUp = goingUpTimes(goingUpTimes>=times.onset(iStim) & goingUpTimes<=times.offset(iStim));
            tmpDown = goingDownTimes(goingDownTimes>=times.onset(iStim) & goingDownTimes<=times.offset(iStim));
            times.frameTimes{iStim} = sort([tmpUp(:); tmpDown(:)], 'ascend');
        end
    end
    
    
    expt.stimTimes(iexp) = times;
    %     expt.stimTimes(ie) = mean(times.onset - Timeline.mpepUDPTimes(startIdx));
    
    expt.phd{iexp} = phd;
    
    
    
    
end

end


% function [ ind ] = conFunStart(str)
%
% ind = contains(str,'StimStart');
%
% end
%
% function [ ind ] = conFunEnd(str)
%
% ind = contains(str,'StimEnd');
%
% end

