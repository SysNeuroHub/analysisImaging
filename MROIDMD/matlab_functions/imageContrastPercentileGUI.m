function [minPct, maxPct] = imageContrastPercentileGUI(img)

    if ndims(img) ~= 2
        error('Input must be a 2D image matrix.');
    end

    img = double(img);
    imgVec = img(:);

    % Initial percentiles
    minPct = 97;
    maxPct = 100;

    % Create figure
    hFig = figure('Name','Percentile Contrast Adjuster',...
                  'NumberTitle','off',...
                  'MenuBar','none',...
                  'ToolBar','none',...
                  'CloseRequestFcn',@closeGUI);

    % Axes
    hAx = axes('Parent',hFig,...
               'Position',[0.1 0.35 0.8 0.6]);

    imagesc(img,'Parent',hAx);
    % colormap(hAx,'gray');
    axis(hAx,'image');

    %% UI

    uicontrol('Style','text','Units','normalized',...
        'Position',[0.02 0.20 0.07 0.05],...
        'String','Min %');

    uicontrol('Style','text','Units','normalized',...
        'Position',[0.02 0.10 0.07 0.05],...
        'String','Max %');

    hMinText = uicontrol('Style','text','Units','normalized',...
        'Position',[0.92 0.20 0.07 0.05],...
        'String',num2str(minPct));

    hMaxText = uicontrol('Style','text','Units','normalized',...
        'Position',[0.92 0.10 0.07 0.05],...
        'String',num2str(maxPct));

    hMinSlider = uicontrol('Style','slider',...
        'Units','normalized',...
        'Position',[0.1 0.20 0.8 0.05],...
        'Min',0,'Max',100,...
        'Value',minPct,...
        'Callback',@updateImage);

    hMaxSlider = uicontrol('Style','slider',...
        'Units','normalized',...
        'Position',[0.1 0.10 0.8 0.05],...
        'Min',0,'Max',100,...
        'Value',maxPct,...
        'Callback',@updateImage);

    handles.hAx = hAx;
    handles.hMinSlider = hMinSlider;
    handles.hMaxSlider = hMaxSlider;
    handles.hMinText = hMinText;
    handles.hMaxText = hMaxText;
    handles.imgVec = imgVec;
    guidata(hFig,handles);

    updateImage();

    uiwait(hFig);

    vals = getappdata(0,'imageContrastPct_values');
    minPct = vals(1);
    maxPct = vals(2);
    rmappdata(0,'imageContrastPct_values');


    %% CALLBACKS

    function updateImage(~,~)

        handles = guidata(hFig);

        minPct = get(handles.hMinSlider,'Value');
        maxPct = get(handles.hMaxSlider,'Value');

        if minPct >= maxPct
            return;
        end

        set(handles.hMinText,'String',sprintf('%.1f',minPct));
        set(handles.hMaxText,'String',sprintf('%.1f',maxPct));

        % Convert percentile → intensity
        climLow  = prctile(handles.imgVec, minPct);
        climHigh = prctile(handles.imgVec, maxPct);

        set(handles.hAx,'CLim',[climLow climHigh]);
    end

    function closeGUI(~,~)
        handles = guidata(hFig);

        minPct = get(handles.hMinSlider,'Value');
        maxPct = get(handles.hMaxSlider,'Value');

        setappdata(0,'imageContrastPct_values',[minPct maxPct]);

        uiresume(hFig);
        delete(hFig);
    end

end