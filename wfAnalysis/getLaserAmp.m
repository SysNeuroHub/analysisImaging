function [laser_interp, tltime] = getLaserAmp(Timeline)
% [laser_interp, tltime] = getLaserAmp(Timeline)
% returns laser amplitude after omitting times for camera exposure


useStimScreen = 1; %use copy of output directory from ignite_wDelay_allOpt

lsr_idx = strcmp({Timeline.hw.inputs.name}, 'laserIn');
laser_raw = Timeline.rawDAQData(:,lsr_idx)'; %[V]
tltime = Timeline.rawDAQTimestamps';

if useStimScreen
    stimScreen_idx =  strcmp({Timeline.hw.inputs.name}, 'stimScreen');
    stimScreen = (Timeline.rawDAQData(:,stimScreen_idx)' > 2.5);
    stimScreenEvent = trace2Event(stimScreen, tltime);
    margin = median(diff(tltime));%0.0005;%[s] %hack
    laserOnTrace = event2Trace(tltime, [stimScreenEvent(:,1)+margin stimScreenEvent(:,2)-margin]);
    laserOnTidx = find(laserOnTrace==1);
else %assume ignite_wDelay_allOpt.ino and parameters hardcoded there
    [strobeOnTimes, strobeOffTimes, strobeDurs] = getStrobeTimes(Timeline, 'alloptrig');
    %frameRate = 60;
    exposureTime = 4; %[ms]
    nRowsOri = 900;
    lineTime =  12.136; %[us]
    bufferTime = 100; %[us]
    sparkDelayTime = 1e-6*(1e3*exposureTime + nRowsOri*lineTime + bufferTime);   %[s]
    
    margin = 0.0005;%[s] %hack
    laserOffTrace = event2Trace(tltime, [strobeOnTimes'-margin strobeOnTimes'+sparkDelayTime+margin]);
    laserOnTidx = find(laserOffTrace==0);
end

laser_interp = interp1(tltime(laserOnTidx), laser_raw(laserOnTidx), tltime); %after omiting camera ON period
