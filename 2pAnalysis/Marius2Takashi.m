function [ImageInfos, roiID] = Marius2Takashi(dat, onlyCells)
%ImageInfos = Marius2Takashi(dat)
%converts result of suite2p to the .ext format compatible with
%Mouse_ExtractFluorescenceChangeNeuropil.m
%
%By default, the function returns only cells, excluding non-cell ROIs
%
%ImageInfos = Marius2Takashi(dat,0)
%returns all ROIs
%
%  This function requires results of new_main(manual curation)
%
%to run this function, make sure to add:
%addpath(genpath('C:\Users\Public\Suite2P')) % add local path to the toolbox
%
%2019-5-26 DS created
%2019-5-31 DS added 2nd input

%OUTPUT:
%ImageInfos
%ImageData
%filterMatrix(_neuropil)
%centerPosition

%TODO:
%DONE check meaning of dat.Fcell is F or dF/F or else
%DONE should multiply Fvalue by npix?
%DONE suppress message from g_getSIstate
%DONE omit non-cell ROIs from output
% convert Fcell from single to double?

if nargin < 2
    onlyCells = 1; %use only ROIs judged as cells
end

try
    ops = dat.ops;
    stat = dat.stat;
catch err
    disp('Input is apparently before manual curation. Use _proc.mat');
end

nTr = length(ops.fsroot{1});

Ly_temp = size(dat.mimg,1); %NG ops.mimg
Lx_temp = size(dat.mimg,2); %NG ops.mimg

if onlyCells
    roiID = find([stat.iscell]); % use only ROIs judged as cells
else
    roiID = 1:length(stat); %use all ROIs
end
stat = stat(roiID); 
nCell = length(roiID);%length(stat);

filterMatrix_temp = zeros(nCell,Ly_temp,Lx_temp);
filterMatrix_neuropil_temp = zeros(nCell,Ly_temp,Lx_temp);
centerPosition_temp = nan(nCell,2);
for iCell = 1:nCell
    filterMatrix_temp(iCell, stat(iCell).ipix) = 1;
    filterMatrix_neuropil_temp(iCell, stat(iCell).ipix_neuropil) = 1;
    centerPosition_temp(iCell,1) = stat(iCell).med(2); %(x)
    centerPosition_temp(iCell,2) = stat(iCell).med(1); %(y)
end

%% translate to image before registration
Ly = size(ops.mimg,1); 
Lx = size(ops.mimg,2); 
filterMatrix = zeros(nCell,Ly,Lx);
filterMatrix(:,ops.yrange,ops.xrange) = filterMatrix_temp;
filterMatrix_neuropil = zeros(nCell,Ly,Lx);
filterMatrix_neuropil(:,ops.yrange,ops.xrange) = filterMatrix_neuropil_temp;
centerPosition(:,1) = centerPosition_temp(:,1) + min(ops.xrange) - 1;
centerPosition(:,2) = centerPosition_temp(:,2) + min(ops.yrange) - 1;

%% extract trace of each trial from concatenated trace (taken from reg2P.m)
Fvalue_conc = dat.Fcell{1}(roiID,:); %OK??
Fvalue_neuropil_conc = dat.FcellNeu{1}(roiID,:);  %OK??
%sp_conc = dat.sp{1}(roiID,:);

fs = ops.fsroot;
nFr = zeros(nTr,1);
startFrame = zeros(nTr,1);
endFrame = zeros(nTr,1);
trialNo = zeros(nTr,1);
iTr = 1; %nbytes = 0;
for k = 1:length(fs)
    % in case different experiments have different numbers of channels
    if ~isempty(ops.expts)
        if ismember(ops.expts(k), getOr(ops, 'expred', []))
            nchannels_expt = ops.nchannels_red;
        else
            nchannels_expt = ops.nchannels;
        end
    else
        nchannels_expt = ops.nchannels;
    end
    if nchannels_expt > 1
        red_mean_expt  = red_mean;
    else
        red_mean_expt  = 0;
    end
    
    for j = 1:length(fs{k})
        nFr(iTr) = nFramesTiff(fs{k}(j).name);
        if isfield(ops,'omitLastFrames') && (ops.omitLastFrames>0) %&& nFr(iTr)>0
            nFr(iTr) = nFr(iTr) - ops.omitLastFrames;
        end
        if mod(nFr(iTr), nchannels_expt) ~= 0
            fprintf('  WARNING: number of frames in tiff (%d) is NOT a multiple of number of channels!\n', j);
        end
        
        
        if iTr == 1
            startFrame(iTr) = 1;
            endFrame(iTr) = nFr(iTr);
        elseif iTr > 1
            startFrame(iTr) = endFrame(iTr-1) + 1;
            endFrame(iTr) = startFrame(iTr) + nFr(iTr) - 1;
        end
        
        [~, filename] = fileparts(fs{k}(j).name);
        filename_parts = textscan(filename,'%s','delimiter','_');
        trialNo(iTr) = str2num(filename_parts{1}{end});
        
        
        state = g_getSIstate(fs{k}(j).name, 1);
        
        ImageData(iTr,1).Fvalue = Fvalue_conc(:,startFrame(iTr):endFrame(iTr));
        ImageData(iTr,1).Fvalue_neuropil = Fvalue_neuropil_conc(:,startFrame(iTr):endFrame(iTr));
        ImageData(iTr,1).trialNo = trialNo(iTr);
        ImageData(iTr,1).relativeTime = calculateRelativeTime(centerPosition, state);
        ImageData(iTr,1).sliceNoList = [stat.iplane]'; %OK??
        ImageData(iTr,1).validityList = int8([stat.iscell]' .*1);%OK??
        ImageData(iTr,1).state = state;
        
        ImageData(iTr,1).DS = ops.DS(startFrame(iTr):endFrame(iTr)); %offsets computed in XY. 27/5/19
        %ImageData(iTr).sp = sp_conc(:,startFrame(iTr):endFrame(iTr));
        
        iTr = iTr + 1;
    end
end
ImageInfos.ImageData = ImageData;
ImageInfos.filterMatrix = filterMatrix;
ImageInfos.centerPosition = centerPosition;
ImageInfos.filterMatrix_neuropil = filterMatrix_neuropil;






