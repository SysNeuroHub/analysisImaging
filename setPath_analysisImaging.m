<<<<<<< HEAD
addpath(genpath('C:\Documents\git\visbox')); %addpath(genpath('C:\Users\Experiment\Documents\MATLAB\visbox'));
addpath(genpath('C:\Documents\git\analysisImaging')); %%addpath(genpath('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\master\Analysis\wfAnalysis\widefield'));    
addpath(genpath('C:\Documents\git\npy-matlab')); %addpath(genpath('C:\npy-matlab'));
addpath(genpath('C:\Documents\git\dsbox'));%addpath(genpath('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\dsbox')); %screen2png
=======
if ispc
    baseDir = 'C:\Users\dshi0006\git';
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
>>>>>>> origin/main
