function [t_concord, input_points, base_points, fliplrCam1, fliplrCam2, ...
    flipudCam1, flipudCam2,cam1_registered, resultFig, scale_recovered, theta_recovered] = ...
    rotate_flip_cams(inputImage, baseImage, method, cmap)


%% make information for rotating inputImage image by visual inspection
% t_concord: rotation info, used for imtransform
% method: transform type (cp2form). default = 'nonreflective similarity'.
%           use 'projective' when 2 exps from 1 animal
%
% input: inputImage, baseImage: 2D or 3D stack of images for registration.
% If 3D, apply the same transformation for each slice in the stack

% should return camera images before thresholding for image registration??

% 2014-5-4 DS added 3rd input
% 2014-5-31 DS added 4th and 5th inputs (default: gray, false)
% 2014-5-10 DS added 10th 11th outputs. allow input and base_images to be 3D
% 2015-12-22 DS added gui sliders to control image min and max
% 2016-03-07 DS added option to pick up landmark points by mouse click
% using ginput

% TODO:
% put location of chosen landmarks on image
% resize multiple images at one time

doOptimize = 0;

if nargin == 4
    if isempty(method)
        method = 'nonreflective similarity';
    end
elseif nargin == 3
    cmap = 'gray';
    if isempty(method)
        method = 'nonreflective similarity';
    end
elseif nargin == 2
    cmap = 'gray';
    method = 'nonreflective similarity';
end
%close;

numInputImage = size(inputImage, 3);
numBaseImage = size(baseImage, 3);

inputCache = linspace(0,1,numInputImage+2);%7/3/16
inputYRange = [inputCache(1:end-1); inputCache(2:end)];
baseCache = linspace(0,1,numBaseImage+2);%7/3/16
baseYRange = [baseCache(1:end-1); baseCache(2:end)];

inputImage(isnan(inputImage)) = 0;
baseImage(isnan(baseImage)) = 0;

for iii = 1:numInputImage
    cacheImage = inputImage(:,:,iii);
    inputImage(:,:,iii) = (inputImage(:,:,iii) - nanmean(cacheImage(:))) ./ nanstd(cacheImage(:));
end
for iii = 1:numBaseImage
    cacheImage = baseImage(:,:,iii);
    baseImage(:,:,iii) = (baseImage(:,:,iii) - nanmean(cacheImage(:))) ./ nanstd(cacheImage(:));
end


inputImage(inputImage>3) = 3;
inputImage(inputImage<-3) = -3;
baseImage(baseImage>3) = 3;
baseImage(baseImage<-3) = -3;

scrsz = get(0,'ScreenSize');
h0 = figure(1);
set(h0,'Position',scrsz);

for iii = 1:numInputImage
    h1(iii) = uipanel('Parent',h0,'position',[0 inputYRange(1,iii) .5 1/numInputImage],...
        'title','inputImage','backgroundcolor','w');
    ax1(iii) = axes('parent',h1(iii));
    imagesc(inputImage(:,:,iii));
    axis equal tight;%grid on;
    pos = get(ax1(iii),'position');
    colorbar('ycolor','k');
    set(ax1(iii),'position',[pos(1) pos(2) pos(3) pos(4)],'xcolor','k','ycolor','k');
    colormap(gca, cmap);grid on;
    
    inputImage_c = squeeze(inputImage(:,:,iii));
    uicontrol(h1(iii),'Style', 'Slider',...
        'Min',min(inputImage_c(:)),'Max',max(inputImage_c(:)),'Value',min(inputImage_c(:)),...
        'Position', [10 20 80 20],...
        'Callback', {@axisMin, ax1(iii)});
    uicontrol(h1(iii),'Style','text',...
        'Position',[20 45 50 20],...
        'String','min')
    
    uicontrol(h1(iii),'Style', 'Slider',...
        'Min',min(inputImage_c(:)),'Max',max(inputImage_c(:)),'Value',max(inputImage_c(:)),...
        'Position', [10 60 80 20],...
        'Callback', {@axisMax, ax1(iii)});
    uicontrol(h1(iii),'Style','text',...
        'Position',[20 85 50 20],...
        'String','max')
    
end

for iii = 1:numBaseImage
    h2(iii) = uipanel('Parent',h0,'position',[0.5 baseYRange(1,iii) .5 1/numBaseImage],...
        'title','baseImage','backgroundcolor','w','visible','on');
    ax2(iii) = axes('parent',h2(iii));
    imagesc(baseImage(:,:,iii));
    axis equal tight;grid minor;
    pos = get(ax2(iii),'position');
    colorbar('ycolor','k');
    set(ax2(iii),'position',[pos(1) pos(2) pos(3) pos(4)],'xcolor','k','ycolor','k');
    colormap(gca,cmap);grid minor;
    
    baseImage_c = squeeze(baseImage(:,:,iii));
    uicontrol(h2(iii),'Style', 'Slider',...
        'Min',min(baseImage_c(:)),'Max',max(baseImage_c(:)),'Value',min(baseImage_c(:)),...
        'Position', [10 20 80 20],...
        'Callback', {@axisMin, ax2(iii)});
    uicontrol(h2(iii),'Style','text',...
        'Position',[20 45 50 20],...
        'String','min')
    
    uicontrol(h2(iii),'Style', 'Slider',...
        'Min',min(baseImage_c(:)),'Max',max(baseImage_c(:)),'Value',max(baseImage_c(:)),...
        'Position', [10 60 80 20],...
        'Callback', {@axisMax, ax2(iii)});
    uicontrol(h2(iii),'Style','text',...
        'Position',[20 85 50 20],...
        'String','max')
end
%close(h0)

validanswer = false;
while ~validanswer
    fliplrCam1 = input('flip left-right inputImage?  ({y}es or {n}o) >>','s');
    for iii = 1:numInputImage
        [inputImage(:,:,iii), validanswer] = flipImage(h0, ax1(iii), inputImage(:,:,iii), fliplrCam1, 'lr',cmap);
        set(ax1(iii),'xcolor','k','ycolor','k');
    end
end

validanswer = false;
while ~validanswer
    flipudCam1 = input('flip up-down inputImage?  ({y}es or {n}o) >>','s');
    for iii = 1:numInputImage
        [inputImage(:,:,iii), validanswer] = flipImage(h0, ax1(iii), inputImage(:,:,iii), flipudCam1, 'ud',cmap);
        set(ax1(iii),'xcolor','k','ycolor','k');
    end
end

validanswer = false;
while ~validanswer
    fliplrCam2 = input('flip lr baseImage?  ({y}es or {n}o) >>','s');
    for iii = 1:numBaseImage
        [baseImage(:,:,iii), validanswer] = flipImage(h0, ax2(iii), baseImage(:,:,iii), fliplrCam2, 'lr',cmap);
        set(ax2(iii),'xcolor','k','ycolor','k');
    end
end

validanswer = false;
while ~validanswer
    flipudCam2 = input('flip ud baseImage?  ({y}es or {n}o) >>','s');
    for iii = 1:numBaseImage
        [baseImage(:,:,iii), validanswer] = flipImage(h0, ax2(iii), baseImage(:,:,iii), flipudCam2, 'ud',cmap);
        set(ax2(iii),'xcolor','k','ycolor','k');
    end
end

validanswer = false;
while ~validanswer
    answer = input('link axis between panels?  ({y}es or {n}o) >>','s');
    if strcmp(answer,'y')
        linkaxes([ax1(:) ax2(:)], 'xy');
        validanswer = true;
    elseif strcmp(answer,'n')
        validanswer = true;
    end
end
zoom on;


validanswer = false;
while ~validanswer
    doRotate = input('rotate inputImage?  ({y}es or {n}o) >>','s');
    
    if strcmp(doRotate, 'y')
        
        %Specifying Control Points
        %h=cpselect(edge(inputImage,'canny'), edge(baseImage,'canny'));
        input_points = []; base_points = [];
        
        %result = 0;%'';
        while 1%result == 0 %~strcmp(result, 'y')
            
            hi = uicontrol(h1(iii),'Style','edit',...
                'Position',[100 20 200 70],...
                'Max',3,'Min',1);%chsnge edit to static string
            set(hi,'String',num2str(input_points)); %put result of previous iteration
            uicontrol(h1(iii),'Style','text',...
                'Position',[100 105 200 20],...
                'String','input points (x, y):')
            set(hi, 'ButtonDownFcn', @addPoints);
            
            
            hb = uicontrol(h2(iii),'Style','edit',...
                'Position',[100 20 200 70],...
                'Max',3,'Min',1);
            set(hb,'String',num2str(base_points)); %put result of previous iteration
            uicontrol(h2(iii),'Style','text',...
                'Position',[100 105 200 20],...
                'String','base points (x, y):')
            set(hb, 'ButtonDownFcn', @addPoints);
            
            %ddd=0;
            hUpdate = uicontrol(h2(iii),'Style','pushbutton',...
                'Position',[350 20 50 40],...
                'String','Update','Callback',@updateEdit);
            set(hUpdate, 'UserData',0);
            
            %             hDone = uicontrol(h2(iii),'Style','pushbutton',...
            %                 'Position',[420 20 50 40],...
            %                 'String','Done','Callback',@doneEdit);
            %             set(hDone, 'UserData',0);
            
            waitfor(hUpdate, 'UserData');
            
            result = ishandle(h0);
            if ~result
                break
            end
            
            
            input_points = str2num(get(hi, 'String'));
            
            
            base_points = str2num(get(hb, 'String'));
            
            
            %Infer Geometric Transformation
            t_concord = cp2tform(input_points,base_points, method);
            
            
            [input_points_rx, input_points_ry] = tformfwd(t_concord,...
                input_points(:,1),input_points(:,2));
            
            if strcmp(method, 'nonreflective similarity')
                Tinv = t_concord.tdata.Tinv;%is this ok?
                ss = Tinv(2,1);
                sc = Tinv(1,1);
                scale_recovered = sqrt(ss*ss + sc*sc);%<1: fitted map is smaller than original anderman's map
                theta_recovered = atan2(ss,sc)*180/pi;
            end
            
            %Transform  Image
            for iii = 1:numInputImage
                cam1_registered(:,:,iii) = imtransform(inputImage(:,:,iii),...
                    t_concord,'XData',[1 size(baseImage,2)], 'YData',[1 size(baseImage,1)]);
            end
            
            %View Registered Image in Context of Orthophoto
            % figure; imshowpair(baseImage,cam1_registered,'blend')
            if exist('resultFig','var')
                clf(resultFig);
            else
                resultFig = figure(2);
                set(resultFig,'Position',scrsz,'visible','on');
            end
            imcache = cam1_registered(:,:,1);
            
            figure(resultFig);
            for iii = 1:numBaseImage
                subplot(numBaseImage, 3, 1+3*(iii-1)); imagesc(imcache, prctile(imcache(:),[15 90]));
                pos = get(gca,'position');
                colorbar;
                set(gca,'position',[pos(1) pos(2) pos(3) pos(4)]);
                axis equal tight;grid minor;title('inputImage rotated');
                hold on;
                plot(base_points(:,1),base_points(:,2),'r*-');
                plot(input_points_rx,input_points_ry,'g*-');
                ax1_result=gca;
                
                imcache = baseImage(:,:,iii);
                subplot(numBaseImage, 3, 2+3*(iii-1)); imagesc(imcache,prctile(imcache(:),[10 90]));
                pos = get(gca,'position');
                colorbar;
                hold on;
                contour(edge(baseImage(:,:,iii),'canny'),1,'r');
                contour(edge(cam1_registered(:,:,1),'canny'),1,'g');
                set(gca,'position',[pos(1) pos(2) pos(3) pos(4)]);
                set(gca,'fontsize',12);grid minor;
                axis equal tight;grid minor;title('baseImage');
                colormap(gca, cmap);
                ax2_result=gca;
                
                subplot(numBaseImage, 3, 3+3*(iii-1));
                imshowpair(baseImage(:,:,iii), cam1_registered(:,:,1), 'scaling','independent');
                set(gca,'fontsize',12);grid minor;
                axis equal tight;
                ax3_result=gca;
            end
            linkaxes([ax1_result ax2_result ax3_result], 'xy');
            
            
            if strcmp(method, 'nonreflective similarity')
                tname = sprintf('scale:%f\nangle:%f', scale_recovered, theta_recovered);
                title(tname);
            end
            
            %% automatic optimization, to refine the alignment ... does not support projective transformation
            if doOptimize
                method_fitgeotrans = method;
                method_fitgeotrans(ismember(method,' '))=[];
                t_concord_test = fitgeotrans(input_points,base_points, method_fitgeotrans);
                
                [optimizer, metric] = imregconfig('multimodal');
                movingReg = imregister2(cam1_registered(:,:,1), baseImage(:,:,iii), 'similarity',optimizer,metric,...
                    'InitialTransformation',t_concord_test);
                figure;
                subplot(121); %this is typically worse than the manual method
                imshowpair(baseImage(:,:,iii), movingReg);
                
                subplot(122); %cropping image helps to improve alignment
                x_crop = round(0.6*size(baseImage,2)):round(0.9*size(baseImage,2));
                y_crop = round(0.3*size(baseImage,1)):round(0.7*size(baseImage,1));
                [~,~,t_concord2] = imregister2(cam1_registered(y_crop,x_crop,1), baseImage(y_crop,x_crop,iii), 'similarity',optimizer,metric,...
                    'InitialTransformation',t_concord_test);
                %input spatial referencing
                RI = imref2d(size(inputImage),[1 size(inputImage,2)],[1 size(inputImage,1)]);
                %output spatial referencing
                Rout = imref2d(size(inputImage),[1 size(baseImage,2)],[1 size(baseImage,1)]);
                cam1_registered2(:,:,iii) = imwarp(inputImage(:,:,iii),RI, ...
                    t_concord2,'OutputView',Rout);
                
                imshowpair(baseImage, cam1_registered2);
            end
            
            %result = input('Are you happy with this registration result? ({y}es or {n}o)>>','s');
            %result = ishandle(h0);
            %result = get(hDone, 'UserData');
            
        end
        validanswer = true;
        
        
    elseif strcmp(doRotate, 'n')
        display('No rotation of inputImage images');
        t_concord = []; input_points = []; base_points = [];
        validanswer = true;
        
    else
        warning('Invalid input. Answer by ''y'' or ''n''');
        validanswer = false;
    end
end

%close all % do not do this here since fig handle is not returned

end

function axisMin(hObject, eventdata, handles)
currentclim = get(handles,'clim');

user_entry = hObject.Value;
if isnan(user_entry)
    return
elseif isempty(user_entry)
    return
else
    set(handles, 'clim', [user_entry currentclim(2)]);
end
end
function axisMax(hObject, eventdata, handles)
currentclim = get(handles,'clim');

user_entry = hObject.Value;
if isnan(user_entry)
    return
elseif isempty(user_entry)
    return
else
    set(handles, 'clim', [currentclim(1) user_entry ]);
end
end

function [inputImage, validanswer] = flipImage(handle_p, ax, inputImage, doFlip, direction, cmap)
if strcmp(doFlip, 'y')
    if strcmp(direction,'lr')
        inputImage = fliplr(inputImage);
    elseif strcmp(direction, 'ud')
        inputImage = flipud(inputImage);
    end
    %ax = axes('parent',handle);
    set(handle_p,'currentaxes',ax);
    imagesc(inputImage, prctile(inputImage(:),[15 97]));
    axis equal tight;
    pos = get(ax,'position');
    colorbar('ycolor','k');
    set(ax,'position',[pos(1) pos(2) pos(3) pos(4)]);
    colormap(gca, cmap);grid minor;
    validanswer = true;
elseif strcmp(doFlip, 'n')
    validanswer = true;
else
    warning('Invalid input. Answer by ''y'' or ''n''');
    validanswer = false;
end
end

function [xx, yy, hp] = addPoints(source, callbackdata)

%Get current figure and axis parameters
oldvals = get(gcf);
oldhold = ishold(gca);
set(gcf,'Pointer','crosshair','doublebuffer','on');
zoom off
[xx,yy] = ginput(1);

currentPoints = str2num(get(source, 'String'));
addedPoints = round([xx yy],1);
set(source, 'String', num2str([currentPoints; addedPoints]));
zoom on

%plot landmark points on the image
plotPoints = [currentPoints; addedPoints];
hold on
hp = plot(plotPoints(:,1),plotPoints(:,2),'ro');
%need to delete previous hp, so only the current plotPoints are visible
%but how to delete previous hp?
hold off

if ~oldhold, hold off; end
end

function source = updateEdit(source, callbackdata)
set(source, 'UserData',1);
end

% function source = doneEdit(source, callbackdata)
% set(source, 'UserData',1);
% end
