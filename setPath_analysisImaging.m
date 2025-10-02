
if ispc
    baseDir = 'C:\Documents\git';
elseif isunix
    baseDir = '/home/daisuke/Documents/git';
end
addpath(genpath(fullfile(baseDir,'visbox'))); %addpath(genpath('C:\Users\Experiment\Documents\MATLAB\visbox'));
addpath(genpath(fullfile(baseDir,'analysisImaging'))); %%addpath(genpath('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\master\Analysis\wfAnalysis\widefield'));    
addpath(genpath(fullfile(baseDir, 'npy-matlab'))); %addpath(genpath('C:\npy-matlab'));
addpath(genpath(fullfile(baseDir, 'dsbox')));%addpath(genpath('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\dsbox')); %screen2png
%addpath('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\master\Analysis\wfAnalysis');
rmpath(genpath(fullfile(baseDir,'analysisImaging/2pAnalysis')));
rmpath(genpath(fullfile(baseDir,'analysisImaging/wfAnalysis/obs')));

