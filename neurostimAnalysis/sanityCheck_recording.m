%this script is meant to be run locally at MU00122377 (data acquisition PC)

% TODO:
% retrieve json file from nsName?
%nsName = 'CJ224.freqTag.185816';
%YYYYMMDD = '20221003';
%
% make plotting function quicker (if possible)?

addpath(genpath('C:\git\analysisImaging'));
addpath(genpath('C:\git\dsbox'));

pathToJsonFile = 'E:\data\2022\10\27\test.oriXYZc.170412_2022-10-27_17-04-29\Record Node 103\experiment1\recording1';

%oeInfo.jsonFile = 'E:\data\2022\10\03\CJ224.freqTag.185816_2022-10-03_18-58-34\Record Node 103\experiment1\recording1\structure.oebin';
oeInfo.jsonFile = fullfile(pathToJsonFile, 'structure.oebin');

%% oephys digital only 
oeInfo.trCh = 1;
oeInfo.camStrobeCh = 7;
oeInfo.ventilatorCh = 4;
oeInfo.expCh = 5;

%% oephys analog
oeInfo.pdCh = 1;

%% load stimulus data
%load(stimFile,'c');

%% retrieve time stamps
disp('Retrieving timestamps in OpenEphys')
OETimes = getOETimes(oeInfo);
