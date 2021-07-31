% compute hemodynamic correction, save Us and Vs to the server
% also compute dF/F, save Us and Vs to the server
% 1/6/20 created from pipelineGCaMP.m for elife paper

addpath('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\master\Analysis\wfAnalysis');
addpath(genpath('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\master\Analysis\wfAnalysis\widefield'));
addpath('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\Matteobox');%mergefigs
addpath('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\dsbox');%screen2png
addpath(genpath('C:\npy-matlab'));

close all;

mouseName = 'TIGRE2GCaMP6s_318';
thisDate = '2021-02-07'; %[datestr(now,'yyyy-mm-dd')];  
thisSeries = 2;
expNums = [2:10 12 15 17];%[1:4];

nSV = 1000;
hemoFreq = 'auto';%
pixSpace = 6; %3
%6 for 0.5x 2x2 binning
%19 for 0.8x 1x1 binning
%24 for 2x 1x1 binning


%where svd is initially saved
%ops.localSavePath = fullfile(rootDrive, 'data', mouseName, thisDate);

%where svd is eventually uploaded

for iexp = 1:length(expNums)
    expPath = fileparts(dat.expFilePath(mouseName, thisDate, thisSeries, expNums(iexp), 'widefield','master'));
    saveVpath = expPath;
    
    quickHemoCorrect(expPath, saveVpath, nSV, hemoFreq, pixSpace);
    
    
    %% dFF
    [U, V, t, mimg] = quickLoadUVt(expPath, nSV, saveVpath);
    [newU, dV] = dffFromSVD(U, V, mimg); %can produce NAN in U
   

    %< TODO when Ublue is in exp directory, change the saving directory
    writeUVtoNPY(newU, dV, fullfile(fileparts(saveVpath), 'svdSpatialComponents_corr_dFF'), ...
        fullfile(saveVpath, 'svdTemporalComponents_corr_dFF'));
    
    copyfile(fullfile(saveVpath, 'svdTemporalComponents_corr.timestamps.npy'),...
        fullfile(saveVpath, 'svdTemporalComponents_corr_dFF.timestamps.npy'));
    
    %% create and save new mimg for quickLoadUVt??
    
    clear U mimg
end