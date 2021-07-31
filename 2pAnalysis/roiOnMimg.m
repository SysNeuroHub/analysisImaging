function cellmap = roiOnMimg(roidata)
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
% 31/3/20 DS created from mapRoidata


cellmap = sum(roidata2map(roidata),3);
imagesc(roidata.image);colormap(gray);hold on
contour(cellmap,'color',[.9 .4 .4]);
axis equal tight


for icell = 1:roidata.num_rois
    text(roidata.centroids(icell,2),roidata.centroids(icell,1),...
        num2str(roidata.s2p_iscell(icell)),'color','r');
end
