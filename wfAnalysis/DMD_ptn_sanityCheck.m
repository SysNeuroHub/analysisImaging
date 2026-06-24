function [f, NG_min] = DMD_ptn_sanityCheck(xfileName, expt, Timeline, stimSequence, DMDIn_state, DMDOut_state)
% [f, NG_min] = DMD_ptn_sanityCheck(xfileName, expt, Timeline, stimSequence, DMDIn_state, DMDOut_state)
% checks DMD pattern order

NG_min = [];

laserTh = 1; %[V]
expt = grabLaserTimesWF(expt,[],laserTh);

tidx = interp1( ...
    Timeline.rawDAQTimestamps, ...
    1:numel(Timeline.rawDAQTimestamps), ...
    expt.laserTimes.onset, ...
    'nearest');
DMDInIdx = DMDIn_state(tidx);
DMDOutIdx = DMDOut_state(tidx);

laserInch = find(strcmp({Timeline.hw.inputs.name}, 'laserIn'));
laserIn_raw= Timeline.rawDAQData(:,laserInch);
laserIn = laserIn_raw > 1;


f = figure('position',[0 0 1500 500]);

if contains(xfileName, 'stimTTLSwitch')
    % 1. DMD in and out idx v p-file
    subplot(221);
    plot(stimSequence.seq, DMDInIdx,'o');xlabel('stimSequence.seq');ylabel('DMD in');
    squareplot;
    
    subplot(222);
    plot(stimSequence.seq, DMDOutIdx,'o');xlabel('stimSequence.seq');ylabel('DMD out');
    squareplot;
    
    [NG_min] = min(find(stimSequence.seq - DMDInIdx ~= 0));
    if ~isempty(NG_min)
        warning('stimSequence and DMD In Index does NOT match');
    end
    
    [NG_min] = min(find(DMDInIdx - DMDOutIdx ~= 0));
    if ~isempty(NG_min)
        warning('DMD In and Out Index does NOT match');
    end
    % [NG_min] = min(find(stimSequence.seq - DMDOutIdx ~= 0));
    
    if ~isempty(NG_min)
        vbox(expt.laserTimes.onset(NG_min), Timeline.rawDAQTimestamps(end),gca, [1 0 0 .5])
    end
end

% 2. DMD in vs out signals over time
subplot(2,4,5);
plot(DMDIn_state, DMDOut_state,'o');xlabel('DMD state IN');ylabel('DMD state OUT');
squareplot;


subplot(2,4,6:8);
plot(Timeline.rawDAQTimestamps, laserIn.*DMDIn_state);hold on
plot(Timeline.rawDAQTimestamps, laserIn.*DMDOut_state);hold on
axis padded

% laserIn_event=trace2Event(laserIn);
%    vbox(Timeline.rawDAQTimestamps(laserIn_event(:,1)), Timeline.rawDAQTimestamps(laserIn_event(:,2)),gca,[.5 .5 1 .2])

legend('DMD state IN','DMD state OUT');
xlabel('time [s]');
ylabel('DMD pattern ID');
