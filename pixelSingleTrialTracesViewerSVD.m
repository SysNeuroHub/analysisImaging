function pixelSingleTrialTracesViewerSVD(U, V, t, mimg, p, stimTimes, ...
    stimSequence, respWin, thisPixel, thisRepeat, subPre)

% 12/11/20 create from showSingleTrialTraces
% TODO: make a version without SVD for 2p data

%p.nstim
%p.nrepeats
%p.filedurs

if nargin < 11
    subPre = 0;
end
if nargin < 10
    thisRepeat = 0; %0 indicates do not show long traces
end

if nargin < 9 || isempty(thisPixel)
    Ypix = size(mimg,1);
    Xpix = size(mimg,2);
    thisPixel = [round(Ypix/2) round(Xpix/2)];
end

% cax = autoCax(allFrames, thisPixel);

f = figure; 

set(f, 'UserData', [thisPixel thisRepeat]);
set(f, 'KeyPressFcn', @(f,k)tcViewerCallbackKey(f, k, U, V, t, ...
    mimg, p, stimTimes, stimSequence, respWin));

redrawFig(U, V, t, mimg, p, stimTimes, ...
    stimSequence, respWin, thisPixel, thisRepeat, subPre);



function tcViewerCallbackClick(f, keydata, U, V, t, mimg, p, stimTimes, ...
    stimSequence, respWin, subPre)
%happens when the top-left image is clicked

figHand = get(f, 'Parent');
ud = get(figHand, 'UserData');
thisPixel = ud(1:2);
thisRepeat = ud(3);

clickX = keydata.IntersectionPoint(1);
clickY = keydata.IntersectionPoint(2);

switch get(f, 'Tag')
    case 'brainImage'
        thisPixel = round([clickY clickX]);        
end

ud(1:2) = thisPixel;

set(figHand, 'UserData', ud);%[thisTimePoint thisPixel thisCond cax]);
redrawFig(U, V, t, mimg, p, stimTimes, ...
    stimSequence, respWin, thisPixel, thisRepeat, subPre);


function tcViewerCallbackKey(f, keydata,  U, V, t, mimg, p, stimTimes, ...
    stimSequence, respWin, subPre)

if ismember(lower(keydata.Key), {'control', 'alt', 'shift'})
    % this happens on the initial press of these keys, so both the Modifier
    % and the Key are one of {'control', 'alt', 'shift'}
    return;
end

%figHand = get(f, 'Parent');
ud = get(f, 'UserData');
thisPixel = ud(1:2);
thisRepeat = ud(3);

%if isequal(keydata.Modifier, {'alt'})
switch lower(keydata.Key)
    case 'uparrow'
        thisRepeat = min(thisRepeat + 1, p.nrepeats);
    case 'downarrow'
        thisRepeat = max(thisRepeat - 1, 1);
end
ud(3) = thisRepeat;
set(f, 'UserData', ud);%[thisTimePoint thisPixel thisCond cax]);
redrawFig(U, V, t, mimg, p, stimTimes, ...
    stimSequence, respWin, thisPixel, thisRepeat, subPre);
  

function redrawFig(U, V, t, mimg, p, stimTimes, ...
    stimSequence, respWin, thisPixel, thisRepeat, subPre)
%(allFrames, cLabels, timePoints, thisPixel, thisTimePoint, thisCond, cax)

%nConditions = size(allFrames,4);

nCols = ceil(sqrt(p.nstim)+1); %extra row is for mimg and long trace
nRows = ceil(p.nstim/(nCols-1))+1;


clf;

% plot the brain image with a marker where the selected pixel is
thisAx = subplot(nRows,nCols,1); 
q = imagesc(mimg); set(q, 'HitTest', 'off');
hold on;
q = plot(thisPixel(2), thisPixel(1), 'ro'); set(q, 'HitTest', 'off');
hold off;
%caxis(cax);
caxis(prctile(mimg(:),[1 99]));
axis equal tight;
title(sprintf('pixel %d, %d selected', thisPixel(1), thisPixel(2)));
set(thisAx, 'ButtonDownFcn', @(f,k)tcViewerCallbackClick(f, k, U, V, t, ...
    mimg, p, stimTimes, stimSequence, respWin, subPre));
% set(thisAx, 'KeyPressFcn', @(f,k)tcViewerCallbackKey(f, k, U, V, t, ...
%     mimg, p, stimTimes, stimSequence, respWin));

set(thisAx, 'Tag', 'brainImage');
colormap(thisAx, 'gray');

% single trial traces of specified pixel
ROItrace = squeeze(svdFrameReconstruct(U(thisPixel(1),thisPixel(2),:),V));
ylimit = prctile(ROItrace,[1 99.8]);
if diff(ylimit) == 0
    ylimit = [ylimit(1) ylimit(1)+1];
end
showSingleTrialTraces(ROItrace, thisRepeat, t, p, stimTimes, stimSequence, ...
    respWin, ylimit, subPre);


function showSingleTrialTraces(ROItrace, irepeat, t, p, stimTimes, stimSequence, ...
    respWin,ylimit, subPre)
% showSingleTrialTraces(ROItrace, irepeat, t, expt, p, stimTimes, stimSequence, respWin)
%
% 22/6/20 created
% 12/11/20 reordered inputs

if nargin < 8
    ylimit = prctile(ROItrace,[1 99]);
end

[~, winSamps, periEvent] = eventLockedAvg(ROItrace', t, stimTimes.onset, stimSequence.seq, respWin);
rectColor = stimID2color(stimSequence.seq, max(stimSequence.seq));

if subPre %12/8/20
    preIdx = find(winSamps<0);
    periEvent = periEvent - mean(periEvent(:,:,preIdx),3);
        ylimit = prctile(periEvent(:),[0.1 99.9]);
end


nCols = ceil(sqrt(p.nstim)+1); %extra row is for mimg and long trace
nRows = ceil(p.nstim/(nCols-1))+1;


%% panel A: concatenated trace of one repeat
% subplot(n_rows,n_rows,1:n_rows);
if irepeat > 0
    eventsInRepeat = (1:p.nstim) + (irepeat-1)*p.nstim;
    subplot(nRows, nCols, 2:nCols);
    
    for ss = eventsInRepeat %size(stimTimes.onset) 1st repeat
        rectangle('position',[stimTimes.onset(ss) ylimit(1) stimTimes.offset(ss)-stimTimes.onset(ss) diff(ylimit)],...
            'edgecolor','none','facecolor',rectColor(ss,:));
        % < stimTimes.offset is buggy...
        hold on;
    end
    tidx = [min(find(t > stimTimes.onset(eventsInRepeat(1))-0.5))...
        :max(find(t < stimTimes.offset(eventsInRepeat(end))+0.5))];
    plot(t(tidx), ROItrace(tidx),'color',[.1 .1 .7]);
    xlim([t(tidx(1)) t(tidx(end))]);%1st repeat
    ylim(ylimit);
    xlabel('Time from exp start [s]');
    ylabel('ROItrace');
    title(['Repeat: ' num2str(irepeat)]);
end

%% panel B: triggered trace of all repeats
avgPeriEvent = nan(p.nstim, length(winSamps));
for istim = 1:p.nstim
    %theseAxes(istim) = subplot(n_rows+1,n_rows,istim+n_rows);
    theseAxes(istim) = subplot(nRows, nCols, istim + nCols);
    
    theseEvents = find(stimSequence.seq == istim);
    
    if ~isempty(theseEvents)
        rectangle('position',[0 ylimit(1) p.pfiledurs(istim) diff(ylimit)],...
            'edgecolor','none','facecolor',stimID2color(istim, max(stimSequence.seq)));
        hold on;
        plot(winSamps, reshape(periEvent(theseEvents,:,:), length(theseEvents),[]),...
            'color',[.5 .5 .5]);
        
        avgTrace = squeeze(mean(periEvent(theseEvents,:,:)));
        avgPeriEvent(istim,:) = avgTrace;
        plot(winSamps, avgTrace, 'k', 'linewidth',2);
        
        if irepeat>0
            thisEvent = intersect(theseEvents, eventsInRepeat);
            plot(winSamps, reshape(periEvent(thisEvent,:,:), 1,[]),...
                'color',[.1 .1 .7]); %12/11/20
        end
    end
    
    title(stimSequence.labels{istim});
    
    xlim([winSamps(1) winSamps(end)]);
    ylim(ylimit);
end
%linkaxes(theseAxes);

subplot(nRows, nCols, p.nstim + nCols + 1);
for istim = 1:p.nstim
    if sum(isnan(squeeze(avgPeriEvent(istim,:)))) == 0
        plot(winSamps, squeeze(avgPeriEvent(istim,:)), ...
            'color', stimID2color(istim, max(stimSequence.seq)),...
            'linewidth',2);
    hold on;
    end
end
xlim([winSamps(1) winSamps(end)]);
ylim(ylimit);
title('avg across repeats')
xlabel('time since stimulus onset [s]');
% caxis([-cax cax]);
% colormap(thisAx, 'jet');colorbar;
% xlabel('azimuth [deg]'); ylabel('elevation [deg]');
% axis equal tight;
