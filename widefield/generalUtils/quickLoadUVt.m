function [U, V, t, mimg, roi] = quickLoadUVt(expPath, nSV, vPath, varargin)
% [U, V, t, mimg, roi] = quickLoadUVt(expPath, nSV[, params])
% loads SVDs
% if hemodynamic correction is already done, load the corrected V
%
% Inputs
%   expPath: fullpath to where SVD data is saved
%   nSV: number of components to load
%   vPath: fullpath to V of SVD
%
% Outputs
%   U: spatial component of SVD
%   V: temporal component of SVD
%   t: time recorded in Timeline time
%   mimg: mean image across frames
%   roi: ROI image specified in pipelineSVDKT, computed from mimg

% %available option:
% movieSuffix = 'blue', useCorrected = 0/1;
% movieSuffix = 'purple', useCorrected = 0;
% movieSuffix = 'corr_dFF', useCorrected = 1;

if nargin < 3 || isempty(vPath)
    vPath = expPath;
end

movieSuffix = 'amber';%'blue';
useCorrected = 1;
if ~isempty(varargin)
    params = varargin{1};
    if isfield(params, 'movieSuffix')
        movieSuffix = params.movieSuffix;
    end
    if isfield(params, 'useCorrected')
        useCorrected = params.useCorrected;
    end
end


expRoot = fileparts(expPath);

disp(['loading svdSpatialComponents_' movieSuffix '.npy'])
try
    U = readUfromNPY(fullfile(expRoot, ['svdSpatialComponents_' movieSuffix '.npy']), nSV);
catch err
    U = readUfromNPY(fullfile(expPath, ['svdSpatialComponents_' movieSuffix '.npy']), nSV);
end

if strcmp(movieSuffix, 'corr_dFF')
    mimg = zeros(size(U,1),size(U,2));
else
    try
        mimg = readNPY(fullfile(expRoot, ['meanImage_' movieSuffix '.npy']));
    catch
        mimg = readNPY(fullfile(expPath, ['meanImage_' movieSuffix '.npy']));
    end
end

% mask out region outside of ROI with NAN 8/7/20
%roiOutside = isnan(U(:,:,1));
roiOutside = (sum(U,3)==0);
mimg(roiOutside)=nan;
roi = ~roiOutside;

timeStampPath = fullfile(expPath, ['svdTemporalComponents_' movieSuffix '.timestamps.npy']);

corrPath = fullfile(vPath, 'svdTemporalComponents_corr.npy');
if exist(corrPath, 'file') && useCorrected
    disp(['loading ' corrPath]);
    %fprintf(1, 'loading corrected temporal components\n');
    V = readVfromNPY(corrPath, nSV);
    t = readNPY(fullfile(expPath, ['svdTemporalComponents_' movieSuffix '.timestamps.npy']));
    Fs = 1/mean(diff(t));

else
    disp(['loading svdTemporalComponents_' movieSuffix '.npy']);
    V = readVfromNPY(fullfile(expPath, ['svdTemporalComponents_' movieSuffix '.npy']), nSV);
    if exist(timeStampPath, 'file')
        t = readNPY(timeStampPath);
        Fs = 1/mean(diff(t));
    else
        disp('timestamps NOT FOUND');
        t = [];
    end
    %V = detrendAndFilt(V, Fs); %17/6/20 commented out
    
end

if length(t)==size(V,2)+1 % happens if there was an extra blue frame at the end
    t = t(1:end-1);
end