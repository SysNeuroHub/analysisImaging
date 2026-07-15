function [minVal, maxVal] = imageContrastGUI(img)

    if ndims(img) ~= 2
        error('Input must be a 2D image matrix.');
    end

    img = normalize_prctile(double(img));
    dataMin = min(img(:));
    dataMax = max(img(:));

    % Create figure
    hFig = figure('Name','Image Contrast Adjuster',...
                  'NumberTitle','off',...
                  'MenuBar','none',...
                  'ToolBar','none',...
                  'CloseRequestFcn',@closeGUI);

    % Axes
    hAx = axes('Parent',hFig,...
               'Position',[0.1 0.35 0.8 0.6]);

    imagesc(img,'Parent',hAx);
    %colormap(hAx,'gray');
    axis(hAx,'image');
    set(hAx,'CLim',[dataMin dataMax]);

    %% ---- UI ELEMENTS ----

    % Static labels
    uicontrol('Style','text','Units','normalized',...
        'Position',[0.02 0.20 0.07 0.05],...
        'String','Min');

    uicontrol('Style','text','Units','normalized',...
        'Position',[0.02 0.10 0.07 0.05],...
        'String','Max');

    % Value display labels
    hMinText = uicontrol('Style','text','Units','normalized',...
        'Position',[0.92 0.20 0.07 0.05],...
        'String',num2str(dataMin,'%.3g'));

    hMaxText = uicontrol('Style','text','Units','normalized',...
        'Position',[0.92 0.10 0.07 0.05],...
        'String',num2str(dataMax,'%.3g'));

    % Min slider
    hMinSlider = uicontrol('Style','slider',...
        'Units','normalized',...
        'Position',[0.1 0.20 0.8 0.05],...
        'Min',dataMin,...
        'Max',dataMax,...
        'Value',dataMin,...
        'Callback',@updateImage);

    % Max slider
    hMaxSlider = uicontrol('Style','slider',...
        'Units','normalized',...
        'Position',[0.1 0.10 0.8 0.05],...
        'Min',dataMin,...
        'Max',dataMax,...
        'Value',dataMax,...
        'Callback',@updateImage);

    % Store handles
    handles.hAx = hAx;
    handles.hMinSlider = hMinSlider;
    handles.hMaxSlider = hMaxSlider;
    handles.hMinText = hMinText;
    handles.hMaxText = hMaxText;
    guidata(hFig,handles);

    uiwait(hFig);

    vals = getappdata(0,'imageContrastGUI_values');
    minVal = vals(1);
    maxVal = vals(2);
    rmappdata(0,'imageContrastGUI_values');


    %% ---- CALLBACKS ----

    function updateImage(~,~)
        handles = guidata(hFig);

        minVal = get(handles.hMinSlider,'Value');
        maxVal = get(handles.hMaxSlider,'Value');

        % Update numeric display
        set(handles.hMinText,'String',num2str(minVal,'%.3g'));
        set(handles.hMaxText,'String',num2str(maxVal,'%.3g'));

        if minVal >= maxVal
            return;
        end

        set(handles.hAx,'CLim',[minVal maxVal]);
    end

    function closeGUI(~,~)
        handles = guidata(hFig);

        minVal = get(handles.hMinSlider,'Value');
        maxVal = get(handles.hMaxSlider,'Value');

        setappdata(0,'imageContrastGUI_values',[minVal maxVal]);

        uiresume(hFig);
        delete(hFig);
    end
end