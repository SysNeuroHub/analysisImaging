function SaveMovieStack(S, iCond, clim, varargin )
% SaveMovieStack(S, iCond, clim, varargin )
% saves a movie
%
% optional inputs;
% frameRate:
% scaleBar
% bgBlack
% mvName
% colorMap
% colorBar
% cbarName
% contour
% 
% returns a video writer class
% a method of StackSet

% TODO: loss-less compression? (video writer class property)
% should return mvobj??
% should specify .avi??
% should save avi, to be compatible with loadAVI??


%default values
frate = 4;
sbar = 'southeast'; %add scale bar ('northwest' 'northeast' 'southwest' 'southeast')
mm = 5; %margin size for scale bar
bgBlack = true; %turn background color black
mvname = 'test.avi'; %name of output avi file
cmap = jet; %specify colormap
cbar = true; %add colorbar
scaleColor = 'k';

if nargin < 3 || isempty(clim);
    % establish a good range
    MedianImg = squeeze(nanmedian(S.Values(:,:,:,iCond),3));
    m = nanmedian(MedianImg(:));
    StdImg = squeeze(nanstd(S.Values(:,:,2:end-1,iCond),[],3));
    StdImg(isinf(StdImg)) = NaN;
    s = nanmean(StdImg(:));
    clim = [ m-2*s m+2*s ];
end

if ~isempty(varargin)
    
    vidx = [];
    for vv = 1 : length(varargin)
        if isstr(varargin{vv})
            vidx = [vidx vv];
        end
    end
    
    %information for frame rate of the movie
    for vv = vidx
        if any(strfind(varargin{vv}, 'frameRate'))
            frate = varargin{vv+1};
            break
        end
    end
    
    %add scale bar
    for vv = vidx
        if any(strfind(varargin{vv}, 'scaleBar'))
            sbar = varargin{vv+1};
             break
        end
    end
    
    %turn background black
    for vv = vidx
        if any(strfind(varargin{vv}, 'bgBlack'))
            bgBlack = varargin{vv+1};
            break
        end
    end
    
    %name of output avi file
    for vv = vidx
        if any(strfind(varargin{vv}, 'mvName'))
            mvname = varargin{vv+1};
            break
        end
    end
    
    %colormap
    for vv = vidx
        if any(strfind(varargin{vv}, 'colorMap'))
            cmap = varargin{vv+1};
            break
        end
    end
    
    %show colorbar
    for vv = vidx
        if any(strfind(varargin{vv}, 'colorBar'))
            cbar = varargin{vv+1};
            break
        end
    end
    
    %add name to colorbar
    for vv = vidx
        if any(strfind(varargin{vv}, 'cbarName'))
            cbarName = varargin{vv+1};
            break
        end
    end
    
    %% add contour map. added on 8/6/14
    for vv = vidx
        if any(strfind(varargin{vv}, 'contour'))
            contourData = varargin{vv+1};
            break
        end
    end
    
end


try
    opengl('software')
    mvobj = VideoWriter(mvname,'motion jpeg avi');
    % Set and view the frame rate.
    mvobj.FrameRate = frate;
    open(mvobj);
    
    
    scrsz = get(0,'ScreenSize');
    ratio = 0.5;
    set(gcf,'Position',[100 100 ratio*scrsz(3) ratio*scrsz(4)],'visible','off');
    
    if bgBlack
        set(gcf,'Color',[0 0 0]);%not complete. also turn the colors of text
        whitebg('k');
    elseif ~bgBlack
        set(gcf,'Color',[1 1 1]);
        whitebg('w'); %added on 8/9/15
    end
    
    colormap(cmap);
    
    for tt = 1:S.nFrames
        imagesc(S.Values(:,:,tt,iCond));
        hold on;
        
        
        if strcmp(sbar, 'northwest')
            line([mm mm+1/S.PixelSize], [mm mm],'color',scaleColor,'linewidth',3);
        elseif  strcmp(sbar, 'northeast')
            line([S.nCols-mm-1/S.PixelSize S.nCols-mm], [mm mm],'color',scaleColor,'linewidth',3);
        elseif strcmp(sbar, 'southwest')
            line([mm mm+1/S.PixelSize], [S.nRows-mm S.nRows-mm],'color',scaleColor,'linewidth',3);
        elseif strcmp(sbar, 'southeast')
            line([S.nCols-mm-1/S.PixelSize S.nCols-mm], [S.nRows-mm S.nRows-mm],'color',scaleColor,'linewidth',3);
        end

        if exist('contourData', 'var')
            contour(contourData,1,'color','b');%should generalize..
        end

        axis equal tight off
        caxis(clim);
        
        if cbar
            pos = get(gca,'position');
            h = colorbar;
            set(gca,'position',[pos(1) pos(2) pos(3) pos(4)]);
            
            if exist('cbarName', 'var')
                ylabel(h, cbarName);
            end
        end
        
        
        tname = sprintf('%.3f [s]',S.TimeVec(tt));
        title(tname,'fontsize',20);
        
        currframe = getframe(gcf);
        writeVideo(mvobj, currframe);%can allow 4d input .. much faster?
        clf;
    end
    
    close(mvobj);
    
catch err
    close(mvobj);
    rethrow(err);
end

close all

end