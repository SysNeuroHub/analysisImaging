function linkcaxes(ax, clim)
%linkcaxes(ax)
% aligns color range across ax
%
%linkcaxes(ax, clim)
% aligns color range across ax with clim
%
%linkcaxes(ax, 'mirror')
% aligns color range across ax with 0 being the center of the colormap

if nargin < 2  ||  strcmp(clim,'mirror')
    for iax = 1:numel(ax)
        clims(:,iax) = ax(iax).CLim;
    end
    clim = [min(clims(:)) max(clims(:))];

    if strcmp(clim,'mirror')
        absvalue = max(abs(clim));
        clim = [-absvalue absvalue];
    end

end


set(ax(:),'clim',clim);