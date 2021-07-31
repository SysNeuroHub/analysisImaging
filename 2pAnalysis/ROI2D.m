function [ROI] = ROI2D(dat, useCells)
%ROIs = ROI2D(dat, useCells)
%returns ROI as a 2D image in the same resolution and size as the original
%
% 16/7/20 DS created from Marius2Craig in Cicada

try
    ops = dat.ops;
    stat = dat.stat;
catch err
    disp('Input is apparently before manual curation (new_main.m). Use _proc.mat');
end

Ly_temp = size(dat.mimg,1); %NG ops.mimg
Lx_temp = size(dat.mimg,2); %NG ops.mimg

if useCells
    roiID = find([stat.iscell]); % use only ROIs judged as cells
else
    roiID = 1:length(stat); %use all ROIs
end
stat = stat(roiID);
nCell = length(roiID);%length(stat);


%% ROI as a 2D mask 
filterMatrix_c = zeros(nCell,Ly_temp,Lx_temp);
for iCell = 1:nCell
    filterMatrix_c(iCell,stat(iCell).ipix) = 1;
end

%% translate to absolute position before registration 17/2/20
filterMatrix = zeros(nCell, dat.ops.Ly, dat.ops.Lx);
filterMatrix(:,dat.ops.yrange, dat.ops.xrange) = filterMatrix_c;

if isfield(ops, 'scaleXY')
    filterMatrix_c =[];%zeros(nCell, newLy, newLx);
    for iCell = 1:nCell
        filterMatrix_c(iCell,:,:) = imresize(squeeze(filterMatrix(iCell,:,:)), 100/ops.scaleXY, 'nearest');
    end
    filterMatrix = filterMatrix_c;
end

if ops.pockelsLineBlank > 0 %26/3/20
    croppedX = ceil(0.5*ops.Lx_ori * ops.pockelsLineBlank/100)+1:floor(ops.Lx_ori - 0.5*ops.Lx_ori * ops.pockelsLineBlank/100);
    filterMatrix_c = zeros(nCell, ops.Ly_ori, ops.Lx_ori);
    filterMatrix_c(:,:,croppedX) = filterMatrix;
    filterMatrix = filterMatrix_c;
end

ROI = shiftdim(filterMatrix,1);