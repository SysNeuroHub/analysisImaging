%NG function [I, ichosen] = redraw_meanimg_test(ax,dat,prop)
function [I] = redraw_meanimg_test(ax,dat,prop)
% used variables in h.dat
% mimg
% dat.map(deleted)
% xlim
% ylim
% axis (renamed from axes2)
% cl.Lx, cl.Ly
% res.iclust

% if h.dat.procmap
%     I = h.dat.mimg_proc(:,:,h.dat.map);
% else    
    I = dat.mimg;%(:,:,dat.map);
% end
 
mu = median(I(:));
sd1 = mean(abs(I(I<mu+1e-7) - mu));
sd2 = 1e-7 + mean(abs(I(I>mu-1e-7) - mu));

%axes(h.axis); 
p=imagesc(ax,I, mu + 5*[-sd1 sd2]);
% imagesc(ax,I);
% caxis(ax, mu + 5*[-sd1 sd2]);
xlim(ax,[dat.xlim]); ylim(ax,[dat.ylim]);
%axis off
axis(ax,'off','equal','tight');
p.ButtonDownFcn = @(o, e) clickHandler(o, e, ax, dat, prop);

colormap(ax,'gray')
%ichosen = dat.ichosen %NG
drawnow
end

function clickHandler(o, e, ax, dat, prop)
x = round(e.IntersectionPoint(1));
y  = round(e.IntersectionPoint(2));

if x>=1 && y>=1 && x<=dat.cl.Lx && y<=dat.cl.Ly && dat.res.iclust(y,x)>0
    dat.ichosen = dat.res.iclust(y, x);
    
    redraw_meanimg_test(ax, dat, prop); %clickHandler will be defined recursively...
    
    str = sprintf('cell #%d\n property %f', dat.ichosen, prop(dat.ichosen));
    title(ax, str);
end
end