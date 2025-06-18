setPath_analysisImaging;

%% experiment
expt.subject = 'yamatotakeru';
expt.expDate = '2025-04-21_1';
expt.expNum = 13;
bklightCtrl = 0;

thisDate = expt.expDate(1:10);
thisSeries = str2num(expt.expDate(12:end));
expPath = fileparts(dat.expFilePath(expt.subject, thisDate, thisSeries, expt.expNum, 'widefield','master'));
saveVpath = expPath;
mpepDir = dat.reposPath('main', 'master');

load(dat.expFilePath(expt.subject, thisDate, thisSeries, expt.expNum, 'Timeline', 'master'));

p = ProtocolLoad_wf(expt.subject,expt.expDate,expt.expNum); %3/6/20


%stim number recorded in p-file
stimSequence = getStimSequence(expt.subject, expt.expDate, expt.expNum);


%stim number according to TTL to DMD
dmdOn = (Timeline.rawDAQData(:,8) > 3);
dmdOn_c = cumsum(dmdOn);
stimSequence_TTL = mod(dmdOn_c, 4);

a(1)=subplot(211);plot(Timeline.rawDAQTimestamps, Timeline.rawDAQData(:,8))
ylabel('DMDin [V]');
a(2)=subplot(212);plot(Timeline.rawDAQTimestamps, stimSequence_TTL);ylim([-.5 4.5]);
ylabel('decoded stim number');xlabel('time[s]');
linkaxes(a(:),'x')
xlim([0 150]); %first
xlim([1460 1520]); %last
