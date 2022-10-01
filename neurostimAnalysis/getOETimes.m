function OETimes = getOETimes(oeInfo, c, makeFig)
%OETimes = getOETimes(oeInfo, c)
%
%INPUT: oeInfo
%.trCh: 1
%.camStrobeCh: 7
%.ventilatorCh: 4
%.expCh: 5
%.pdCh: 1
%.jsonFile: (fullpath to the oephys json file)

%OUTPUT:
%OETimes.expOn(Off)Times
%OETimes.stimOn(Off)TimesDigi
%OETimes.stimOn(Off)TimesPD
%OETimes.ventOn(Off)Times
%OETimes.camOn(Off)Times
%OETimes.stimOn(Off)Times: stimOnset times after adjustment

if nargin < 3
    makeFig = 1;
end

nrRepeats = c.nrTrials;%    c.blocks.nrRepeats;

[OETimes.expOnTimes, OETimes.expOffTimes] = getDigitalTimesOE(oeInfo.jsonFile, oeInfo.expCh);
[OETimes.stimOnTimesDigi, OETimes.stimOffTimesDigi] = getDigitalTimesOE(oeInfo.jsonFile, oeInfo.trCh);
[OETimes.stimOnTimesPD, OETimes.stimOffTimesPD, pd_ds, pdOn, t_ori, t_ds] = processContinuousDiodeOE(oeInfo.jsonFile, oeInfo.pdCh);
[OETimes.ventOnTimes, OETimes.ventOffTimes] = getDigitalTimesOE(oeInfo.jsonFile, oeInfo.ventilatorCh);
[camOnTimes_c1, camOnTimes_c2] = getDigitalTimesOE(oeInfo.jsonFile, oeInfo.camStrobeCh);
OETimes.camOnTimes = sort([camOnTimes_c1;camOnTimes_c2]);



if isempty(OETimes.stimOffTimesPD)
    OETimes.stimOffTimesPD = OETimes.camOnTimes(end); %HACK
end
if isempty(OETimes.stimOffTimesDigi)
    OETimes.stimOffTimesDigi = OETimes.camOnTimes(end); %HACK
end


if makeFig
    ax(1)=subplot(411);
    try
        vbox(OETimes.stimOnTimesDigi, OETimes.stimOffTimesDigi);
    end
    vline([OETimes.expOnTimes OETimes.expOffTimes]);
    ylabel('stimOnTimes Digital');
    title('getOETimes');
    
    ax(2)=subplot(412);
    plot(t_ds,pd_ds./prctile(pd_ds,99))
    try
        vbox(OETimes.stimOnTimesPD, OETimes.stimOffTimesPD);
    end
    vline([OETimes.expOnTimes OETimes.expOffTimes]);
    ylabel('stimOnTimes PD');
    
    ax(3)=subplot(413);
    vbox(OETimes.ventOnTimes, OETimes.ventOffTimes);
    vline([OETimes.expOnTimes OETimes.expOffTimes]);
    ylabel('ventilation');
    
    ax(4)=subplot(414);
    vbox(camOnTimes_c1, camOnTimes_c2);
    vline([OETimes.expOnTimes OETimes.expOffTimes]);
    ylabel('camera triggering');
    
    linkaxes(ax(:));
end

if ~isempty(OETimes.stimOnTimesPD)
    [OETimes.stimOnTimes, OETimes.stimOffTimes] = checkAdjustTrTimes(OETimes.stimOnTimesPD, ...
        OETimes.stimOffTimesPD, OETimes.expOnTimes, nrRepeats, makeFig);
else
    [OETimes.stimOnTimes, OETimes.stimOffTimes] = checkAdjustTrTimes(OETimes.stimOnTimesDigi, ...
    OETimes.stimOffTimesDigi, OETimes.expOnTimes, nrRepeats, makeFig);
end 
