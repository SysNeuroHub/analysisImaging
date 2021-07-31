%NG function [I, dat] = redraw_figure_test(ax, dat, prop) 
function [I] = redraw_figure_test(ax, dat, prop)
% show an HSV image, values of which are determined by prop
% this function is called after selecting new cell (housed in
% h.dat.F.ichosen)
%
% prop must be [0 1]

%fields in h.dat
% cl.Ly
% cl.Lx
% cl.rands (will be replaced with some functional property)
% cl.vmap (deleted)
% xlim
% ylim
% axis (renamed from axes2)
% ichosen (renamed from F.ichosen)
% stat
%   ipix
%   lambda (deleted)
%   lam
%   iscell (deleted)
%

%TODO return ichosen
if nargin < 3
    prop = dat.cl.rands;
end

Sat1     =  ones(dat.cl.Ly, dat.cl.Lx);
H1              = zeros(dat.cl.Ly, dat.cl.Lx);

[~, V1] = ...
    getviclust(dat.stat, dat.cl.Ly,  dat.cl.Lx, dat.ichosen);

iclust = dat.res.iclust;

%iselect     = iclust1==dat.ichosen; %pixels of clicked ROI
iselect     = iclust ==dat.ichosen; %pixels of clicked ROI
Sat1(iselect)= 0;
% H1(iclust1>0)   = dat.cl.rands(iclust1(iclust1>0));
%H1(iclust1>0)   = prop(iclust1(iclust1>0));
H1(iclust>0)   = prop(iclust(iclust>0));%outside ROI as black

I = hsv2rgb(cat(3, H1, Sat1, V1));
I = min(I, 1);
p = imagesc(ax, I);
xlim(ax,[dat.xlim]); ylim(ax,[dat.ylim]);
axis(ax, 'equal','tight','off');
p.ButtonDownFcn = @(o, e) clickHandler(o, e, ax, dat, prop);

% ichosen = dat.ichosen; %correct?
end
 
function [iclust1, V1] = getviclust(stat, Ly, Lx, ichosen)
% vmap = 'var' or 'unit'
% ichosen: cell id of selected by mouse clicking
% OUTPUT
%   iclust1:
%   V1: selected cell will be 1 (white)

iclust1 = zeros(Ly, Lx);
V1      = zeros(Ly, Lx);

for j = 1:numel(stat)
    ipix    = stat(j).ipix;
    lambda   = stat(j).lam;
    
    if ichosen==j
        inew = true(numel(ipix), 1);
    else
        inew    = lambda(:)>(V1(ipix) + 1e-6);
    end
    
    L0      = stat(j).lam(inew);
    V1(ipix(inew))      = L0;
    iclust1(ipix(inew)) = j;
 
end
mV = mean(V1(V1>0));
V1 = V1/mV;
end

function clickHandler(o, e, ax, dat, prop)
x = round(e.IntersectionPoint(1));
y  = round(e.IntersectionPoint(2));

if x>=1 && y>=1 && x<=dat.cl.Lx && y<=dat.cl.Ly && dat.res.iclust(y,x)>0
    dat.ichosen = dat.res.iclust(y, x);
    
    redraw_figure_test(ax, dat, prop); %clickHandler will be defined recursively...
    
    str = sprintf('cell #%d\n property %f',dat.ichosen, prop(dat.ichosen));
    title(ax, str);
end
end