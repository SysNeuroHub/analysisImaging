function [ expt ] = grabSparseFrames( expt )
%[ expt ] = grabSparseFrames( expt )
% retrieve time of stimulus presentation from p-file
% OUTPUT:
%  expt.stimFrames(stimNum).ss
%  expt.stimFrames(stimNum).ny
%  expt.stimFrames(stimNum).nx

% make sure to add:
% addpath('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\master\visbox\Stimulus');
% m file for the stimulus protocol
% addpath('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\master\stimFiles\current');

%acceptable xfiles:
% stimSparseNoiseUncorrAsync


%% commented out 10/7/20
% timelineFile = expt.timelineFile;
% 
% % basicVS_autoLFR(timelineFile{1});
% 
% load(timelineFile);%{1}
% 
% %expNum = expt.expNum; %expr
% 
% 
% i=1;
% while i < Timeline.mpepUDPCount
%     if strfind(Timeline.mpepUDPEvents{i,1}, 'BlockStart')
%         infoSaveString = Timeline.mpepUDPEvents{i,1};
%         expInfo = textscan(infoSaveString(12:end), '%s', 'delimiter',' ');
%         
%         %SetDefaultDirs;
%         %Protocol = ProtocolLoad(cell2mat(expInfo{1}(1)),str2num(cell2mat(expInfo{1}(2))),str2num(cell2mat(expInfo{1}(3))));
%         
         Protocol = ProtocolLoad_wf(expt.subject,expt.expDate,expt.expNum); %3/6/20
% 
%         i=Timeline.mpepUDPCount;
%     else
%         i = i+1;
%     end
% end


% fsHardware = dir(sprintf('%s/*hardwareInfo.mat', expt.timelineDir));
% load(sprintf('%s/%s', expt.timelineDir{1}, fsHardware(1).name));
thisDate = expt.expDate(1:10);
thisSeries = str2num(expt.expDate(12:end));
load(dat.expFilePath(expt.subject, thisDate, thisSeries, expt.expNum, 'hw-info','master'));

myScreenInfo.windowPtr = NaN; % so we can call the stimulus generation and it won't try to display anything
expt.stimFrames = [];
for stimNum = 1:Protocol.nstim
    ss = eval([Protocol.xfile(1:end-2) '(myScreenInfo, Protocol.pars(:,stimNum));']);
    stim_screen = cat(3,ss.ImageTextures{:});
    ny = size(stim_screen,1);
    nx = size(stim_screen,2);
    
    expt.stimFrames(stimNum,1).ss = ss;
    expt.stimFrames(stimNum,1).ny = ny;
    expt.stimFrames(stimNum,1).nx = nx;
end

end

