function [laser_interp, tltime] = getLaserAmp(Timeline)
% [laser_interp, tltime] = getLaserAmp(Timeline)
% returns laser amplitude after omitting times for camera exposure
% currently assuming ignite_wDelay_allOpt with fixed parameters

lsr_idx = strcmp({Timeline.hw.inputs.name}, 'laserIn');
laser_raw = Timeline.rawDAQData(:,lsr_idx)'; %[V]
tltime = Timeline.rawDAQTimestamps';
[strobeOnTimes, strobeOffTimes, strobeDurs] = getStrobeTimes(Timeline, 'alloptrig');

%ignite_wDelay_allOpt.ino
%frameRate = 60;
exposureTime = 4; %[ms]
nRowsOri = 900;
lineTime =  12.136; %[us]
bufferTime = 100; %[us]
sparkDelayTime = 1e-6*(1e3*exposureTime + nRowsOri*lineTime + bufferTime);   %[s]

margin = 0.0005;%[s] %hack
laserOffTrace = event2Trace(tltime, [strobeOnTimes'-margin strobeOnTimes'+sparkDelayTime+margin]);
laserOnTidx = find(laserOffTrace==0);

laser_interp = interp1(tltime(laserOnTidx), laser_raw(laserOnTidx), tltime); %after omiting camera ON period
