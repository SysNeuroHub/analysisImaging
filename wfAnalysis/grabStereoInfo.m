function expt = grabStereoInfo(expt)
%expt = grabStereoInfo(expt)
%output:
%expt.stereoInfo created by saveStereoInfo

thisDate = expt.expDate(1:10);
thisSeries = str2num(expt.expDate(12:end));
expPath = fileparts(dat.expFilePath(expt.subject, thisDate, thisSeries, expt.expNum, 'widefield','master'));
saveVpath = expPath;

filename_stereo = sprintf('%s_%d_%s_stereoInfo.mat',expt.expDate, expt.expNum, expt.subject);
fullpath_stereo = fullfile(expPath, filename_stereo);

if exist(fullpath_stereo, 'file')
    load(fullpath_stereo, 'stereoInfo');
else
    disp('StereoInfo not found. Make from scratch');
    stereoInfo = saveStereoInfo(expt);
end
expt.stereoInfo = stereoInfo;