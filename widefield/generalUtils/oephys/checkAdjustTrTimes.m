function [trOnTimes, trOffTimes] = checkAdjustTrTimes(trOnTimes, trOffTimes, expOnTimes, nrRepeats, makeFig)
%[stimOnTimesPD, stimOffTimesPD] = checkAdjustTrTimes(stimOnTimesPD, stimOffTimesPD, expOnTimes, nrRepeats)

if nargin < 5
    makeFig = true;
end

if makeFig
    figure;
    ax(1)=subplot(211);
    %vbox(trOnTimes, trOffTimes)
    vline(trOnTimes, ax(1), '-',[1 0 0]);hold on;
    vline(trOffTimes, ax(1), '-',[0 0 1]);
    vline(expOnTimes);
    xlabel('time');
    ylabel('before adjustment');
    title('checkAdjustTrTimes');
end

%% sanity check trial onset timestamps. Adjust if necessary
%first trial onset appears after the experiment start?
if isempty(expOnTimes)
    error('experiment onset test: experiment onset NOT detected');
elseif numel(expOnTimes)>1
    error('experiment onset test: Multiple experiment onset detected');
elseif expOnTimes > trOnTimes(1)
    disp('Stimulus onset detected earlier than experiment onset');
    disp('removed the first trial onset&offset');
    trOnTimes = trOnTimes(2:end);
    trOffTimes = trOffTimes(2:end);
else
    disp('PASSED experiment onset test:  Stimulus onset came after experiment onset');
end

%#trial onsets and offsets are equal?
if numel(trOnTimes) ~= numel(trOffTimes)
    [trOnTimes,trOffTimes] = adjustOnOffTimes(trOnTimes,trOffTimes);
    disp('adjusted trial on/off numbers');
%     disp(['trial on/off numbers test: # trial onsets detected (' num2str(numel(trOnTimes)) ') larger than # trial offsets detected (' num2str(numel(trOffTimes)) ')'])
%     disp('removed the last trial onset');
%     trOnTimes = trOnTimes(1:end-1);
% elseif numel(trOnTimes) < numel(trOffTimes)
%     disp(['trial on/off numbers test: # trial onsets detected (' num2str(numel(trOnTimes)) ') smaller than # trial offsets detected (' num2str(numel(trOffTimes)) ')'])
%     if trOnTimes(1)>=trOffTimes(1)
%         trOffTimes = trOffTimes(2:end);
%         disp('removed the first trial offset');
%     else
%         error('maybe initial trial onset NOT detected? dont know what to do');
%     end
else
    disp(['PASSED trial on/off numbers test: '  num2str(numel(trOnTimes)) 'trial onsets detected']);
end
    
%recorded #trial equal to one specified in CIC?
if numel(trOnTimes) > nrRepeats
    error('trial numbers in CIC test: Too many stimulus onsets detected. Consider adjusting threshold');
elseif numel(trOnTimes) < nrRepeats
    warning(['trial numbers in CIC test: ' num2str(nrRepeats-numel(trOnTimes)) 'missed trials.']);
else 
    disp('PASSED trial numbers in CIC test: #detected trials matches #specified trials');
end

if makeFig
    ax(2)=subplot(212);
    vbox(trOnTimes, trOffTimes)
    vline(expOnTimes);
    xlabel('time');
    ylabel('after adjustment');
    linkaxes(ax(:));
end
end

function [xstart,xend] = adjustOnOffTimes(xstart, xend)
%used when xstart and xend have unequal number of events
%copied from vbox.m

if numel(xstart)<numel(xend)
    if xstart(1) > xend(1)
        xend = xend(2:end);
    else
        xend = xend(1:end-1);
    end
elseif numel(xstart)>numel(xend)
    if xstart(end) > xend(end)
        xstart = xstart(1:end-1);
    else
        xstart = xstart(2:end);
    end
end
end