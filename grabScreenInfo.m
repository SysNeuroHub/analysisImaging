function expt = grabScreenInfo(expt, SIDir)
%expt = grabScreenInfo(expt, SIDir)
% retrieves screenInfo of the specified experiment, put it as expt.screenInfo
%12/10/20 created from grabStimTimesWF.m

if nargin < 2 || isempty(SIDir)
    SIDir = '\\storage.erc.monash.edu.au\shares\MNHS-dshi0006\Subjects';
end

filename_SI = sprintf('%s_%d_%s_hardwareInfo.mat',expt.expDate, expt.expNum, expt.subject);
screenInfo_filename = fullfile(SIDir, expt.subject, expt.expDate, num2str(expt.expNum), filename_SI);
load(screenInfo_filename, 'myScreenInfo');

expt.screenInfo = myScreenInfo;