function dat = smooth_movie(dat, ops)
% dat = smooth_movie(dat, ops)
% applies smoothing only for estimating dispersion in image registration
% using ops.smooth_time_space
% if length(smooth_time_space) == 1
%   apply smoothing in time. smoth_time_space defines std [frames]
% if length(smooth_time_space) == 2
%   apply smoothing in time then space, with the same smoothing factor between x & y
% if length(smooth_time_space) == 3
%   apply smoothing in time and space, with different smoothing factor in x & y

% 2020 Oct copied from Suite2p/registration
% called in registration_offsets.m

smooth = reshape(ops.smooth_time_space, 1, []);
switch length(smooth)
    case 1
        smDims = 3;
        smSigs = smooth;
        if smooth <= 0
            smooth = [];
        end
    case 2
        smSigs = smooth([2 2 1]);
        smDims = find(smSigs > 0);
        smSigs = smSigs(smDims);
        if isempty(smDims)
            smooth = [];
        end
    case 3
        smSigs = smooth([2 3 1]);
        smDims = find(smSigs > 0);
        smSigs = smSigs(smDims);
        if isempty(smDims)
            smooth = [];
        end
end
if ~isempty(smooth)
    dat = my_conv2(double(dat), smSigs, smDims);
    dat = int16(dat); %commented out 20/10/20
end