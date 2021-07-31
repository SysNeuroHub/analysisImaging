function stereoInfo = saveStereoInfo(expt, movieSuffix, refImage)
%stereoInfo = saveStereoInfo(expt, movieSuffix, refImage)
% 15/10/20 created from buildStereotaxicInfo_test (shimaoka elife 2019)

%TODO:
%save directory one above  exp
%copy existing stereoInfo

if nargin < 3
    refImage = []; %to be implemented
end
if nargin < 2 || isempty(movieSuffix)
    movieSuffix = 'blue';
end

%INPUT:
% expt.subject
% expt.expDate
% expt.expNum
%
% movieSuffix: 'blue','purple'

%OUTPUT:
% stereoInfo:
% bregma [y x]
% lambda [y x]
% yaw [deg]
% pitch [deg]
% roll [deg]

thisDate = expt.expDate(1:10);
thisSeries = str2num(expt.expDate(12:end));
expPath = fileparts(dat.expFilePath(expt.subject, thisDate, thisSeries, expt.expNum, 'widefield','master'));
saveVpath = expPath;

filename_stereo = sprintf('%s_%d_%s_stereoInfo.mat',expt.expDate, expt.expNum, expt.subject);
fullpath_stereo = fullfile(expPath, filename_stereo);

params.movieSuffix = movieSuffix;
params.useCorrected = 0;
[~,~,~, mimg] = quickLoadUVt(expPath, 1, saveVpath, params);


stereoInfo = [];

%% find bregma position in pixels
[bregmax, bregmay] = findPixLocation(mimg, 'bregma');
[lambdax, lambday] = findPixLocation(mimg, 'lambda');
stereoInfo.bregma = [bregmay bregmax];
stereoInfo.lambda = [lambday lambdax];

%% find yaw angle
stereoInfo.yaw = showYawRotateImage(mimg);

%% input pitch and roll angle
prompt = {'pitch angle [deg]','roll angle [deg]'};
name = 'Input experiment parameters';
numlines =1;
defaultanswer = {'0','20'};
options.Resize = 'on';
options.WindowStyle = 'normal';
options.Interpreter = 'tex';
answer = inputdlg(prompt, name, numlines, defaultanswer, options);
stereoInfo.pitch = str2num(answer{1});
stereoInfo.roll = str2num(answer{2});


%% upload to server
save(fullpath_stereo,'stereoInfo');

close all




%% inline functions

    function outputyawAngle = showYawRotateImage(inputImage)
        
       
            outputyawAngle = 0;
        inputImage_ori = inputImage;
        %         thMin = prctile(double(inputImage(:)), 1);
        %         thMax = prctile(double(inputImage(:)), 99);
        
        % Create a figure and axes
        f = figure('Visible','off');
        drawImages(inputImage);
        
        % Create slider
        sld = uicontrol('Style', 'slider',...
            'Min',-90,'Max',90,'Value',0,...
            'sliderstep',[0.5/180 10/180],...
            'Position', [400 20 120 20],...
            'Callback', @setAngle);
        
        % Add a text uicontrol to label the slider.
        txt = uicontrol('Style','text',...
            'Position',[400 45 120 20],...
            'String','yaw angle');
        
        % Make figure visble after adding all components
        f.Visible = 'on';
        % This code uses dot notation to set properties.
        % Dot notation runs in R2014b and later.
        % For R2014a and earlier: set(f,'Visible','on');
        
        uiwait(f);%wait until a figure is closed
        
        function setAngle(source,callbackdata)
            outputyawAngle = source.Value;
            
            inputImage = imrotate(inputImage_ori, outputyawAngle);
            
            drawImages(inputImage)
        end
        
    end

    function drawImages(inputImage)
        
        
        imagesc(inputImage);
        
        hold on
        %contour(maskIdx2D,[.5 .5],'r');
        %contour(edge(inputImage(:,:,1),'canny'),'color','r');
        %contour(edge(inputImage(:,:,2),'canny'),'color','g');
        hold off
        
        axis equal tight
        grid on
    end


    function [bregmaxpix, bregmaypix] = findPixLocation(inputImage, thistext)
        % Create a figure and axes
        f = figure('Visible','off');
        %ax = axes;
        
        
%         inputImage = imadjust(inputImage, [0 1]); %scale images from 0 to 1
        crange = prctile(double(inputImage(:)), [1 99]);
        
        
       
        bregmaxpix = 1;
        bregmaypix = 1;
        imhandle = drawWithPlaid(inputImage, crange, bregmaxpix, bregmaypix);
        
        
        
        % Create slider to modulate max/min of background image
        sld2 = uicontrol('Style', 'slider',...
            'Min',min(inputImage(:)),'Max',max(inputImage(:)),'Value',prctile(double(inputImage(:)),1),...
            'Position', [10 100 100 30],...
            'Callback', @setMin);
        
        sld3 = uicontrol('Style', 'slider',...
            'Min',min(inputImage(:)),'Max',max(inputImage(:)),'Value',prctile(double(inputImage(:)),99),...
            'Position', [10 130 100 30],...
            'Callback', @setMax);
        
        set(f, 'WindowButtonDownFcn', @addPosition);
        
        % Make the uibuttongroup visible after creating child objects.
        f.Visible = 'on';
        
        msg = ['double click position of ' thistext ', then close the figure'];
        h = msgbox(msg);
        
        uiwait(f);%wait until a figure is closed
        
        
        function imhandle = drawWithPlaid(inputImage, crange, bregmaxpix, bregmaypix)
            
            %inputImage = imadjust(inputImage, crange, [0 1]); %for RGB image
            
            imhandle = imagesc(inputImage);
            
            hold on
            %             contour(edge(inputImage(:,:,1),'canny'),'color','r');
            %             contour(edge(inputImage(:,:,2),'canny'),'color','g');
            
            plot(bregmaxpix, bregmaypix, 'mo');
%             xgrid = 1/Exps.Cam.MmPerPixel * linspace(-5,5,11)+bregmaxpix;
%             ygrid = 1/Exps.Cam.MmPerPixel * linspace(-5,5,11)+bregmaypix;
%             
%             set(gca,'xtick',xgrid,'xticklabel',linspace(-5,5,11),...
%                 'ytick',ygrid,'yticklabel',-linspace(-5,5,11));
            grid on;
            
            caxis(crange);
            axis equal tight
            hold off
        end
        
        
        function setMin(source,callbackdata)
            crange(1) = source.Value;
            
            imhandle = drawWithPlaid(inputImage, crange, bregmaxpix, bregmaypix);
        end
        
        function setMax(source,callbackdata)
            crange(2) = source.Value;
            
            imhandle = drawWithPlaid(inputImage, crange, bregmaxpix, bregmaypix);
        end
        
        function addPosition(source, callbackdata)
            %[myobj,xs,ys] = freehanddraw_s(imhandle, callbackdata, gca);
            [bregmaxpix, bregmaypix] = ginput(1);
            
            imhandle = drawWithPlaid(inputImage, crange, bregmaxpix, bregmaypix);
        end
        
    end


end


