

% expInfo.nsName = 'CJ229.oriXYZc.170024';%'CJ229.runPassiveMovies.024114';
% expInfo.expID = 5; %21;
% expInfo.subject = 'CJ229';
% expInfo.date = '202210310';%'20221101';

expInfo.subject = 'CJ224';
expInfo.date = '20221004';
expInfo.nsName = 'CJ224.runPassiveMovies.033059';
expInfo.expID = 19;

procParam.rescaleFac = 0.25;
procParam.cutoffFreq = 0.01;%0.1;
procParam.lpFreq = 2; %1

rebuildImageData = false;
makeMask = false;
uploadResult = true;

%downsample of stim & imaging data
dsRate = 1;%[Hz]

moviePath = 'Z:\Shared\Daisuke\natural\nishimoto2011';

% gabor bank filter 
gaborBankParamIdx.cparamIdx = 1;
gaborBankParamIdx.gparamIdx = 2;
gaborBankParamIdx.nlparamIdx = 1;
gaborBankParamIdx.dsparamIdx = 1;
gaborBankParamIdx.nrmparamIdx = 1;



%% save image and processed data
saveDir = fullfile(oeOriDir, fullOEName);
imProc_fullpath = fullfile(saveDir, imDataName);
    
if exist(imProc_fullpath,'file') % && 
    load(imProc_fullpath,'imageProc');
else
    imageProc = saveImageProcess(expInfo, procParam, rebuildImageData,...
    makeMask, uploadResult);
end

observed = prepareObserved(imageProc, dsRate);
% save as timetable?
% TT = timetable(seconds(TimeVec_ds_c), observed);
% writetimetable(TT, timeTableSaveName);


%% prepare stimulus (~5h in my PC)
S_fin = saveStimData(moviePath, imageProc.cic, imageProc.stimInfo, dsRate, ...
    gaborBankParamIdx, 0);

dirPref = getpref('nsAnalysis','dirPref');
%oeOriServer = dirPref.oeOriServer; %direct ethernet connection
%saveDir = fullfile(oeOriDir, fullOEName);

%% save gabor filter output as .mat
stimSaveName = fullfile(dirPref.saveDirBase, expDate,...
    ['stimData_' regexprep(expDate, filesep, '_') '_' num2str(expInfo.expID) '.mat']);
save( stimSaveName, 'TimeVec_stim_cat', 'S_fin');%, 'energyModelParams');


