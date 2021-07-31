function cellmap = mapRoidata(roidata, showFig)
% cellmap=mapRoidata(roidata)
% shows superimposed image of ROIs and its ID in suite2p notation
% also returns the image data
%
% INPUT: "roidata" format 
%  roidata.rois_ibar
%  roidata.s2p_iscell
%  roidata.centroids
% 
% see also: load_S2P_rois
%
% 30/5/2019 DS

if nargin < 2
    showFig = 1;
end

dims = size(roidata.image);
cellmap = sum(roidata2map(roidata),3);

if showFig
    imagesc(cellmap);axis equal tight;colormap(1-gray);
    hold on
    for icell = 1:roidata.num_rois
        text(roidata.centroids(icell,2),roidata.centroids(icell,1),...
            num2str(roidata.s2p_iscell(icell)),'color','r');
    end
    grid minor
end