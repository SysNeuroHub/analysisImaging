% compute hemodynamic correction, save Us and Vs to the server
% also compute dF/F, save Us and Vs to the server
% 1/6/20 created from pipelineGCaMP.m for elife paper

% addpath('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\Matteobox');%mergefigs
% addpath(genpath('C:\Users\dshi0006\git\dsbox'));
% addpath(genpath('C:\Users\dshi0006\git\analysisImaging'));
% addpath(genpath('C:\npy-matlab'));
setPath_analysisImaging;

close all;

mouseName = 'apollo';
thisDate = '2025-05-21'; %[datestr(now,'yyyy-mm-dd')];  
thisSeries = 1;
expNums = [2];%[1:4];

nSV = 2000; %NG with 1000 components
hemoFreq = [7 9];% [5 9];%
pixSpace = 1; %6
%6 for 0.5x 2x2 binning
%19 for 0.8x 1x1 binning
%24 for 2x 1x1 binning


%where svd is initially saved
%ops.localSavePath = fullfile(rootDrive, 'data', mouseName, thisDate);

%where svd is eventually uploaded

for iexp = 1:length(expNums)
    expPath = fileparts(dat.expFilePath(mouseName, thisDate, thisSeries, expNums(iexp), 'widefield','master'));
    saveVpath = expPath;
    
    quickHemoCorrect_binTest(expPath, saveVpath, nSV, hemoFreq, pixSpace);
    
    
%     %% dFF
%     [U, V, t, mimg] = quickLoadUVt(expPath, nSV, saveVpath);
%     [newU, dV] = dffFromSVD(U, V, mimg); %can produce NAN in U
%    
% 
%     %< TODO when Ublue is in exp directory, change the saving directory
%     writeUVtoNPY(newU, dV, fullfile(fileparts(saveVpath), 'svdSpatialComponents_corr_dFF'), ...
%         fullfile(saveVpath, 'svdTemporalComponents_corr_dFF'));
%     
%     copyfile(fullfile(saveVpath, 'svdTemporalComponents_corr.timestamps.npy'),...
%         fullfile(saveVpath, 'svdTemporalComponents_corr_dFF.timestamps.npy'));
    
    %% create and save new mimg for quickLoadUVt??
    
    clear U mimg
end