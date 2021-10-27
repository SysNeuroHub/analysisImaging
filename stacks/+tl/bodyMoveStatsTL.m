function [speed_median, speed_std, SPEED, TIME, TIME_c] = ...
    bodyMoveStatsTL(Exps, stimList, repeatList, TLstring, windowSizeSec, drawSummary, timeFromVsSt)
%
% [speed_median, speed_std] = bodyMoveStatsTL(Exps, iStim)
% returns median and std of speed [m/s] (in a trial) of locomotion of every trials
%
% [...] = bodyMoveStatsTL(Exps, iStim, repeatList)
% returns speed of specified trials
%
% [...] = bodyMoveStatsTL(Exps, iStim, repeatList, TLstring)
% allows to specify which Timeline channel to use (default:'rotaryEncoder')
%
%  [...] = bodyMoveStatsTL(Exps, iStim, repeatList, TLstring, windowSize)
% allows to pecify length of smoothing window in seconds
%
% [speed_median, speed_std, SPEED, TIME] = bodyMoveStatsTL(...) also returns cell
% array of speed of each trial, and time from vs onset
% 
% [speed_median, speed_std, SPEED, TIME, TIME_c] = bodyMoveStatsTL(...) also returns cell
% array of time in Timeline time

%SetDefaultDirs;

% 2014-10-21 DS : allow iStim to be a vector (stimlist)
% 2014-11-05 DS : allow repeatList to be a cell, specifying repeatList for
% each stimulus
% 2015-12-08 DS : first 3 inputs are replaced with Exps (Exps.animal, Exps.iseries, Exps.iexp)
% 2017-7-3 DS : added 5th input

%timeFromVsSt = false;%24/4/15

if nargin < 7 % DS on 5/6/15
    timeFromVsSt = true;
end



if nargin < 6
    drawSummary = false;
end
if nargin < 5
    windowSizeSec = 0.05;
end
if nargin < 4 || isempty(TLstring)
    TLstring = 'rotaryEncoder';
end



if nargin < 1 % input parameters via a pop-up window
    prompt={'Animal name:',...
        'Serie num',...
        'Experiment num',...
        'Stimulus num',...
        'Trial num (empty all)'...
        'Name in Timeline',...
        'Temporal smoothing size',...
        'Summary figure? (y/n)',...
        };
    name='Input experiment pars';
    numlines=1;
    defaultanswer={sprintf('M%6.0f_SD',str2num(datestr(date,'yymmdd'))),...
        '1','1','1','','rotaryEncoder','50','y'};
    
    options.Resize='on';
    options.WindowStyle='normal';
    options.Interpreter='tex';
    
    answer = inputdlg(prompt,name,numlines,defaultanswer,options);
    
    animal = answer{1};
    iSeries = str2num(answer{2});
    iExp = str2num(answer{3});
    
    stimList = str2num(answer{4});
    repeatList = str2num(answer{5});
    TLstring = answer{6};
    windowSizeSec = str2num(answer{7});
    if strcmp(str2num(answer{8}), 'y')
        drawSummary = true;
    else
        drawSummary = false;
    end
    
    Exps = [];
    Exps.animal = animal;
    Exps.iseries = iSeries;
    Exps.iexp = iExp;
end

p = ProtocolLoad(Exps);

%4/9/2015 DS for parallel computing
% p = ProtocolLoad(Exps.animal, Exps.iseries, Exps.iexp, 'donotload',...
%         '\\zserver.ioo.ucl.ac.uk\Data\trodes');

ppr = 1024;%number of pulses per revolution %KÜBLER - 05.2400.1122.0100
ballDiameter = 0.197; %[m] diamter of treadmil
mpp = pi*ballDiameter / ppr / 4; %[m/pulse] on treadmil %2014/11/28 fixed

if isempty(stimList)
    stimList = p.nstim;
end

if ~exist('repeatList','var')
    repeatList = [];
end
if isempty(repeatList) || nargin < 5
    repeatList = 1:p.nrepeats;%TO DO: change this for each exp
end

if ~iscell(repeatList)
    for jjj = 1:length(stimList)
        repeatList_cache{jjj} = repeatList;
    end
    repeatList = repeatList_cache;
end

tlname = sprintf('//zserver/Data/expInfo/%s/%d/%d/%d_%d_%s_Timeline.mat',...
    Exps.animal,Exps.iseries, Exps.iexp, Exps.iseries, Exps.iexp, Exps.animal);
disp('bodyMoveStatsTL:Loading Timeline...')
load(tlname);
srate = Timeline.hw.daqSampleRate;%sampling rate of daq [Hz]

for jjj = 1:length(stimList)
    
    iStim = stimList(jjj);
    
    for iii = 1:length(repeatList{jjj})
        
        iTr = repeatList{jjj}(iii);
        
        [stIdx, enIdx] = ...
            tl.extractTimeTL(Timeline, Exps.animal, Exps.iseries, Exps.iexp, iStim, iTr);
        
        if timeFromVsSt
            %vsStTime = Timeline.rawDAQTimestamps(vsStIdx);
            
            vsStTime = tl.visStimStart(Timeline, Exps.animal, Exps.iseries, Exps.iexp, ...
                iStim, iTr); %14/8/15
        else
            vsStTime = 0;
        end
        
        %extract data sequence (name specified by camString)
        [rotarySequence] = tl.extractSequenceTL(Timeline, stIdx, enIdx, TLstring);
        
        %         %correct reverse movement. 1/8/2014...does not work to M141014_SD
        %         difRotarySequence = diff(rotarySequence);
        %         idx = find(difRotarySequence > 4e9);
        %         difRotarySequence(idx) = -1;
        %         idx = find(difRotarySequence < -4e9);
        %         difRotarySequence(idx) = +1;
        %         %rotarySequence = [rotarySequence(1); cumsum(difRotarySequence)];
        %         rotarySequence = cumsum(difRotarySequence);
        
        % correct reverse movement.
        idx = find(rotarySequence>1e8);
        maxValue = max(rotarySequence(idx)); %28/11/14
        rotarySequence(idx) = rotarySequence(idx) - maxValue; %28/11/14
        
        windowSize = round(windowSizeSec * srate);
        
        rotaryFilt = filtfilt(ones(1,windowSize)/windowSize,1,double(rotarySequence));
        speed = srate*diff(rotaryFilt); %[pulses/sec]
        speed = mpp * speed; %[m/pulse]*[pulses/sec]
        
        speed_std(iii) = std(speed);
        speed_median(iii) = median(speed);
        
        if nargout > 2
            SPEED{jjj, iii} = speed;
        end
        if nargout > 3
            %tentative..
            %TIME{jjj, iii} = Timeline.rawDAQTimestamps(stIdx+1:enIdx-1)' - vsStTime;
            TIME{jjj, iii} = Timeline.rawDAQTimestamps(stIdx:enIdx-1)' + 0.5*Timeline.hw.samplingInterval - vsStTime;
        end
        if nargout > 4
            TIME_c{jjj, iii} = Timeline.rawDAQTimestamps(stIdx:enIdx-1)';
        end
    end
end