function roidata = load_S2P_rois(roipos, mimage)
% [rois_ibar,centroids,areas,s2p_iscell] = load_S2P_rois(roipos, dims)
% returns rois_ibar of CELLS (not ROIs)
%
% INPUT:
% roipos(roinumber).xpos (absolute position)
% roipos(roinumber).ypos (absolute position)
% roipos(roinumber).iscell copy of dat.stat.iscell
%
% OUTPUT:
% roidata(roinumber).rois_ibar 
% roidata(roinumber).s2p_iscell
% roidata(roinumber).rois_ibar
% roidata(roinumber).roi_idcs
% roidata(roinumber).num_rois
% roidata(roinumber).roi_status
% roidata(roinumber).image
%
% necessary information for this function
% dat.stat.iscell
% dat.stat.xpix, dat.stat.ypix
% dat.ops.xrange, dat.ops.yrange
%
% 30/5/2019 DS created from registers2p.m

dims = size(mimage);
if isfield(roipos, 'iscell')
    s2p_iscell = find([roipos(:).iscell]);
else
    s2p_iscell = 1:length(roipos);
end
%rois_ibar = zeros(sum(cellfun(@(x) numel(x),{dat.stat(s2p_iscell).ypix}')),4);
rois_ibar = zeros(sum(cellfun(@(x) numel(x),{roipos(s2p_iscell).ypos}')),4);
prev_id = 1;
centroids = zeros(numel(s2p_iscell),2);
areas = centroids(:,1);
for r = 1:numel(s2p_iscell)
    
    ypos = roipos(s2p_iscell(r)).ypos;
    xpos = roipos(s2p_iscell(r)).xpos;

    %min_x = dat.ops.xrange(1);
    %min_y = dat.ops.yrange(1);
    centroids(r,:) = mean([ypos xpos],1);
    
    ypos_validmin = find(ypos >= 1);
    ypos_validmax = find(ypos <= dims(1));%2/4/20
    yidx_valid = intersect(ypos_validmin, ypos_validmax);
    
    xpos_validmin = find(xpos >= 1);
    xpos_validmax = find(xpos <= dims(2));
    xidx_valid = intersect(xpos_validmin, xpos_validmax);
    
    xyidx_valid = intersect(yidx_valid, xidx_valid);
    
    these_idcs = sub2ind(dims, ypos(xyidx_valid), xpos(xyidx_valid));
    %             these_idcs = sub2ind(dims, dat.stat(s2p_iscell(r)).ypix+min_y-1, ...
    %                 dat.stat(s2p_iscell(r)).xpix+min_x-1);
    
    
    areas(r) = numel(these_idcs);
    blank_im = zeros(dims);
    blank_im(these_idcs) = 1;
    b = bwperim(blank_im);
    b = find(b);
    border = 0 * these_idcs;
    border(ismember(these_idcs,b)) = 1;
    area = ~border; %1 + (0 * these_idcs);
    %rois_ibar(prev_id:prev_id+numel(these_idcs)-1,:) = [these_idcs border area  (r + (0 * these_idcs))];
    rois_ibar(prev_id:prev_id+numel(these_idcs)-1,:) = ...
        [these_idcs border area  (s2p_iscell(r) + (0 * these_idcs))]; % Added 20170511 HD to fix ROI indexing
    prev_id = prev_id + numel(these_idcs);
end

rois_ibar = rois_ibar(rois_ibar(:,1)>0,:);%29/5/19

% added on 31/5/2019
roidata.s2p_iscell = s2p_iscell;
roidata.centroids = centroids;
roidata.rois_ibar = rois_ibar;
roidata.roi_idcs = unique(roidata.rois_ibar(:,4),'stable');
roidata.num_rois = numel(unique(roidata.rois_ibar(:,4)));
roidata.roi_status = [ones(size(roidata.rois_ibar,1),1) ...
    zeros(size(roidata.rois_ibar,1),1)];
%roidata.image = mimage;

end


