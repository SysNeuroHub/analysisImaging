function showSingleTrialTraces(ROItrace, irepeat, t, p, stimTimes, stimSequence, respWin,ylimit)
% showSingleTrialTraces(ROItrace, irepeat, t, expt, p, stimTimes, stimSequence, respWin)
%
% 22/6/20 created
% 12/11/20 reordered inputs

if nargin < 8
    ylimit = prctile(ROItrace,[1 99]);
end

[avgPeriEvent, winSamps, periEvent] = eventLockedAvg(ROItrace', t, stimTimes.onset, stimSequence.seq, respWin);
rectColor = stimID2color(stimSequence.seq, max(stimSequence.seq));


n_rows            = ceil(sqrt(p.nstim));

eventsInRepeat = (1:p.nstim) + (irepeat-1)*p.nstim;

%% panel A: concatenated trace of one repeat
subplot(n_rows+1,n_rows,1:n_rows);
for ss = eventsInRepeat %size(stimTimes.onset) 1st repeat
    rectangle('position',[stimTimes.onset(ss) ylimit(1) stimTimes.offset(ss)-stimTimes.onset(ss) diff(ylimit)],...
        'edgecolor','none','facecolor',rectColor(ss,:));
    hold on;
end
tidx = [min(find(t > stimTimes.onset(eventsInRepeat(1))-0.5))...
    :max(find(t < stimTimes.offset(eventsInRepeat(end))+0.5))];
plot(t(tidx), ROItrace(tidx),'k');
xlim([t(tidx(1)) t(tidx(end))]);%1st repeat
ylim(ylimit);
xlabel('Time from exp start [s]');
ylabel('ROItrace');
title(['Repeat: ' num2str(irepeat)]);

%% panel B: triggered trace of all repeats
for istim = 1:p.nstim
    theseAxes(istim) = subplot(n_rows+1,n_rows,istim+n_rows);
    
    theseEvents = find(stimSequence.seq == istim);
        
    rectangle('position',[0 ylimit(1) p.pfiledurs(istim) diff(ylimit)],...
        'edgecolor','none','facecolor',stimID2color(istim, max(stimSequence.seq)));
    hold on;
    plot(winSamps, reshape(periEvent(theseEvents,:,:), length(theseEvents),[]));
    plot(winSamps, squeeze(avgPeriEvent(istim,:,:)), 'k', 'linewidth',2);
    title(stimSequence.labels{istim});
    
    xlim([winSamps(1) winSamps(end)]);
    ylim(ylimit);
end
%linkaxes(theseAxes);
xlabel('time since stimulus onset [s]');