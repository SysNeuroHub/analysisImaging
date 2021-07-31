function map = roidata2map(roidata)
%map = roidata2map(roidata)
% returns 2Dmap of each ROI
%
% INPUT: roidata, result of load_S2P_rois
% OUTPUT: 2D map of ROI (x-y-roi)
%
% 2019-5-31 DS

dims = size(roidata.image);
map = zeros(dims(1),dims(2),roidata.num_rois);
for icell = 1:roidata.num_rois
    map_c = zeros(dims);
    map_c(roidata.rois_ibar(find(roidata.rois_ibar(:,4) == ...
        roidata.s2p_iscell(icell)),1)) = 1;
    map(:,:,icell) = map_c;
end