function drawTopDownCtx(addLabels, addColors, cortexOnly, bregmaXY, contourColor)
% drawTopDownCtx(addLabels, addColors, cortexOnly, bregmaXY, contourColor)
% draw top down view of cortex borders
% requires curve fitting toolbox (smooth)
% inputs:
% addLabels: if true, show labels of each area (default false)
% addColors: if true, show color from CCF (default: false)
% cortexOnly: if true, show only isocortex (default: true_
% bregmaXY: specify location of bregma. By default, use the estimate by NS
% contourColor: color of the contour (default: 'w' meaning white)

% 2018-3-30 DS created from script_TopDownCortex.m

load('dataForTopDownCtx','tdAnn','st');

lwidth = 0.5; %18/10/18

if nargin < 1 || isempty(addLabels)
    addLabels = false;
end
if nargin < 2 || isempty(addColors)
    addColors = false;
end
if nargin < 3 || isempty(cortexOnly)
    cortexOnly = true;
end
if nargin < 4 || isempty(bregmaXY)
    bregma = allenCCFbregma;
    bregmaXY = [bregma(3) bregma(1)];
end
if nargin < 5
    contourColor = 'w';
end

if addColors
    cmap = allen_ccf_colormap('2017'); % this colormap is the official allen-chosen colors
end

%figure;
hold on;
axis image; %axis off; set(gca, 'YDir','reverse');
uAnn = unique(tdAnn(:));

if cortexOnly 
    isoctxID = st.id(strcmp(st.acronym, 'Isocortex'));
    uAnn = uAnn(cellfun(@(x)contains(x, ['/' num2str(isoctxID) '/']), st.structure_id_path(uAnn)));
else
    uAnn = uAnn(uAnn~=1 & uAnn~=1107 & uAnn~=334 & uAnn~=381 & uAnn~=1143); % exclude root, onl, RSPd, MOB, bic
end

se = strel('disk', 3);


if ~isempty(bregmaXY)
    scale = 1e-2; %mm/pix
    xaxis = scale*((0:size(tdAnn,2)+1) - bregmaXY(1));
    yaxis = scale*((0:size(tdAnn,1)+1) - bregmaXY(2));
else
    xaxis = 0:size(tdAnn,2)+1;
    yaxis = 0:size(tdAnn,1)+1;
end
[XAXIS,YAXIS] = meshgrid(xaxis,yaxis);
    
for u = numel(uAnn):-1:1
    if uAnn(u)==380
        indInSlice = tdAnn==uAnn(u) | tdAnn==381 | tdAnn==1107; %special case for olfactory, which otherwise looks funny
    else
        indInSlice = tdAnn==uAnn(u);
    end
    indInSlice = imerode(imdilate(imdilate(imerode(indInSlice, se),se),se),se); 

    % pad with zeros otherwise contours on the edge do a funny thing
    blankCanvas = zeros(size(indInSlice,1)+2, size(indInSlice,2)+2);
    blankCanvas(2:end-1,2:end-1) = indInSlice; indInSlice = blankCanvas;

    [c,~] = contour(XAXIS, YAXIS, double(indInSlice), [0.5 0.5]); %this is slow
    ii = 1;
    while ii<size(c,2)
        n = c(2,ii);
        if n>100
            x = c(1,ii+1:ii+n);
            y = c(2,ii+1:ii+n);
            x = smooth(x,50,'loess');
            y = smooth(y,50,'loess');
            x = [x;x(1)];
            y = [y;y(1)];
            if addColors
                plot(x,y, 'LineWidth', lwidth, 'Color', cmap(uAnn(u),:));
            else
                plot(x,y, 'LineWidth', lwidth, 'Color', contourColor);
            end
            
            % using a fill looks nicer but doesn't save to pdf very well
            %fill(x,y, cmap(uAnn(u),:), 'FaceAlpha', 0.2, 'EdgeAlpha', 0);
            
            if addLabels
                mnx = mean(x);%-20; 
                mny = mean(y);%-8;
                str = st.acronym{uAnn(u)}; 
                if str(end)=='1'; str = str(1:end-1); end
                ah = annotation('textbox', [0.1, 0.1, 0.1, 0.1], 'String', ...
                    str, 'FontSize', 7, 'LineStyle', 'none');
                set(ah, 'Parent', gca); set(ah, 'Position', ...
                    [mnx,mny, 0.1, 0.1]);
            end
            
        end
        ii = ii+n+1;
    end
    drawnow;
end
