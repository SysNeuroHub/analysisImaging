function [stimOnframes, stimOR] = grabORSequenceFrames( expt )
% [stimOnframes, stimOR] = grabORSequenceFrames( expt )
%21/8/20 created from mmn_CR_DS1_2018_11_27_6

thisDate = expt.expDate(1:10);
thisSeries = str2num(expt.expDate(12:end));
load(dat.expFilePath(expt.subject, thisDate, thisSeries, expt.expNum, 'hw-info','master'));

myScreenInfo.windowPtr = NaN; % so we can call the stimulus generation and it won't try to display anything

Protocol = ProtocolLoad_wf(expt.subject,expt.expDate,expt.expNum); %3/6/20

for stimNum = 1:Protocol.nstim
    ss = eval([Protocol.xfile(1:end-2) '(myScreenInfo, Protocol.pars(:,stimNum));']);
    
    stimOnframes{stimNum} = [1 1+find(diff(ss.ImageSequence(1,:),[],2) == -1)];
    %stimOnframes{stimNum} = [1 1+find(diff(ss.ImageSequence(1,:),[],2) ~= 0)];
    stimOR{stimNum} = ss.Orientations(1,stimOnframes{stimNum});
end