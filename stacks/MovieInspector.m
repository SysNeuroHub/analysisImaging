function varargout = MovieInspector(varargin)
% MovieInspector allows you to inspect any movie
% 
% MovieInspector(mmm) allows you to inspect the movie mmm, which is 
% nrows by ncolumns by nframes
%
% MovieInspector(mmm,roi) lets you specify one or more regions of interest.
% Default = [], meaning no ROI.
%
% MovieInspector(mmm,roi, fs) lets you specify the sampling rate in Hz
% (DEFAULT: []).
%
% MovieInspector(mmm,roi, fs, clim) lets you specify the clipping limits
% (DEFAULT: []).
%
% MovieInspector(mmm,roi, fs, clim, TitleString) lets you specify the title
% string
%
% MovieInspector(mmm,roi, fs, clim, TitleString, TimeVec) lets you specify
% time axis
%
% MovieInspector(mmm,roi, fs, clim, TitleString, TimeVec, overlayImage)
% lets you specify transparency image on the movie
%
% ****** Example 1 (reading a monochromatic AVI) ****** 
%
% % Read the AVI
% m = VideoReader('ExampleAVI.avi');
% nframes = m.FrameRate*m.Duration;
% mmm = zeros(m.Height,m.Width,nframes );
% for iframe = 1:nframes
%     mmm(:,:,iframe) = mean(read(m,iframe),3);
% end
% % Crop it
% mmm([ 1:29 376:end ],:,:  ) = []; mmm(:, [ 1:113 471:end ],:) = [];
% 
% % Inspect it
% MovieInspector(mmm)
% 
% ****** Example 2 (reading a monochromatic Tiff stack) ****** 
%
% InfoImage=imfinfo('ExampleTiffStack.tif');
% mImage=InfoImage(1).Width;
% nImage=InfoImage(1).Height;
% nFrames=length(InfoImage);
% TifLink = Tiff('ExampleTiffStack.tif', 'r');
% 
% mmm=zeros(nImage,mImage,nFrames,'uint16');
% for iFrame=1:nFrames
%    TifLink.setDirectory(iFrame);
%    mmm(:,:,iFrame)=TifLink.read();
% end
% TifLink.close();
% 
% mmm(:,:,1) = []; % lose the first frame, it is weird
% 
% % Inspect it
% MovieInspector(double(mmm));
%
% 2013-03 Matteo Carandini
% 2013-04 MC added support for frame rate (3rd argument)
% 2013-10 MC added TitleString
% 2014-01 DS added an option to save video (parameters should be
% 2014-03 MC copied from MovieInspector but put frames as 3rd dimension
% 2015-02 DS added 6th(TimeVec) and 7th(overlayImage) inputs 
% controllable from outside)
% 2015-03 DS added 7th input(colormap)
% 2015-03 DS fixed bug of range to show histogram
% 2015-03 DS fixed bug of color of axes of time trace
% 2015-04 DS debugged 2nd input (roi), added slider for speed control, 
%            merged ">" and "||" buttons, change line properties of roi rectangles 
% 2015-04 DS debugged so it is compatible with matlab 2013b or older
% 2016-01 MC removed bug that made timecourses clickable, changed slider range

% MOVIEINSPECTOR MATLAB code for MovieInspector.fig
%      MOVIEINSPECTOR, by itself, creates a new MOVIEINSPECTOR or raises the existing
%      singleton*.
%
%      H = MOVIEINSPECTOR returns the handle to a new MOVIEINSPECTOR or the handle to
%      the existing singleton*.
%
%      MOVIEINSPECTOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MOVIEINSPECTOR.M with the given input arguments.
%
%      MOVIEINSPECTOR('Property','Value',...) creates a new MOVIEINSPECTOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MovieInspector_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MovieInspector_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MovieInspector


% TO DO
% Swap ">" and ">>" buttons to "ll" during playing
% 
% Show multiple conditions at the same time


% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MovieInspector_OpeningFcn, ...
                   'gui_OutputFcn',  @MovieInspector_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

% End initialization code - DO NOT EDIT


% --- Executes just before MovieInspector is made visible.
function MovieInspector_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MovieInspector (see VARARGIN)

handles.mmm = varargin{1}; %tensor data

if length(varargin)<8 % DS on 2015-2-24
    cmap = [];
else
    cmap = varargin{8};
end

if length(varargin)<7 % DS on 2015-2-24
    overlayImage = [];
else
    overlayImage = varargin{7};
end

if length(varargin)<6 % DS on 2014-7-28
    TimeVec = [];
else
    TimeVec = varargin{6};
end

if length(varargin)<5
    TitleString = 'MovieInspector';
else 
    TitleString = varargin{5};
end

if length(varargin)<4
    clim = [];
else
    clim = varargin{4};
end

if length(varargin)<3
    fs = [];
else
    fs = varargin{3};
end

if length(varargin)<2
	roi = [];
else
    roi = varargin{2};
end

set(handles.axTimeCourse,'NextPlot','add');

handles.nFrames = size(handles.mmm,3);

if ~isempty(TimeVec)
    handles.tt = TimeVec;
elseif ~isempty(fs)
    handles.tt = (1:handles.nFrames)/fs;
else
    handles.tt = (1:handles.nFrames);
end

if isempty(clim)
    % MC introduced ff on 2014-02-08
    ff = isfinite(handles.mmm(:));
    clim = [min(handles.mmm(:)) max(handles.mmm(:))]; % too slow: prctile(handles.mmm(ff),[.5 99.5]);
end

iFrame = 1;

set(handles.figMovieInspector,'Name',TitleString);

handles.TimeMark = plot( handles.axTimeCourse,handles.tt([iFrame iFrame]),clim, 'k-' );
set(handles.axTimeCourse,'xlim',handles.tt([1 end]),'ylim',clim,'box','off','color','none',...
    'Xcolor','k','ycolor','k','fontsize',10);
 
set(handles.axMovie, 'XColor','k','YColor','k');
axes(handles.axMovie);
if isempty(cmap)
    colormap bone
else
    colormap(cmap);
end
if isempty(overlayImage)
    handles.Image = imagesc(squeeze(handles.mmm(:,:,iFrame)),clim);
else
    handles.Image = imagesc(squeeze(handles.mmm(:,:,iFrame)),...
        'alphadatamapping','scaled','alphadata',overlayImage, clim);
end
hold on
% colorbar
axis image

handles.TimePlots = [];
handles.roiPoints = [];
for iroi = 1:size(roi,1)
    roi_cache = imrect (handles.axMovie, roi(iroi,:));
    AddOneROI(hObject,handles,roi_cache);
    handles = guidata(hObject); % refresh the handles
end

set(handles.Image, 'ButtonDownFcn',...
     @(hObject,eventdata)MovieInspector('axMovie_ButtonDownFcn',hObject,eventdata,guidata(hObject)) );
 
handles.CurrentFrame = 1;

% upon closing we should delete this timer (we probably do already)
handles.Timer = timer(...
    'Period',0.1,...
    'ExecutionMode','fixedRate',...
    'UserData',hObject,...
    'TimerFcn',@(hObject,eventdata)MovieInspector('ShowNextFrame',hObject,eventdata));

% Flag whether to save the movie or not
handles.saveMov = 0; %2014-7-28 DS added


% Choose default command line output for MovieInspector
handles.output = hObject;

opengl('software')
mvname = 'test.avi';
handles.mvobj = VideoWriter(mvname,'motion jpeg avi');

% Set parameters for the video.
handles.mvobj.FrameRate = 5;%10;
handles.mvobj.Quality = 50;%100;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MovieInspector wait for user response (see UIRESUME)
% uiwait(handles.figMovieInspector);


% --- Outputs from this function are returned to the command line.
function varargout = MovieInspector_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in pbPlayPause.
function pbPlayPause_Callback(hObject, eventdata, handles)
if strcmp(handles.Timer.Running,'off') %when movie is NOT being played, start
    pbPlay_Callback(hObject, eventdata, handles);
elseif handles.Timer.Running %handles.Playing %when movie is being played, stop
    pbPause_Callback(hObject, eventdata, handles);
end

% --- Executes on button press in pbPlayPause.
function pbPlay_Callback(hObject, eventdata, handles)
% hObject    handle to pbPlayPause (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%handles.figMovieInspector.Children(length(handles.figMovieInspector.Children) - 2).String = '||';
set(handles.pbPlayPause, 'String', '||');
guidata(hObject, handles);
if handles.Timer.Running, stop(handles.Timer); end
%set(handles.Timer,'Period',0.1);
start(handles.Timer);


% --- Executes on button press in pbPlayPause.
function pbPause_Callback(hObject, eventdata, handles)
% hObject    handle to pbPause (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%handles.figMovieInspector.Children(length(handles.figMovieInspector.Children) - 2).String = '>';
set(handles.pbPlayPause, 'String', '>');
guidata(hObject, handles);
stop(handles.Timer);

if handles.saveMov
    handles = guidata(hObject);
    handles.saveMov = 0;
    guidata(hObject, handles);
    
    disp('Recording mode OFF');
    pause(0.5);%to suppress the following error: "OBJ must be open before writing video"
    
    close(handles.mvobj);
end

% --- Executes on button press in pbToStart.
function pbToStart_Callback(hObject, eventdata, handles)
% hObject    handle to pbToStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ShowFrame(hObject, handles,1);

% Utility function to show a frame
function ShowFrame(hObject, handles, iFrame)
set(handles.TimeMark,'xdata',handles.tt([iFrame iFrame]));
set(handles.Image,'CData', squeeze(handles.mmm(:,:,iFrame)) );  
handles.CurrentFrame = iFrame;
% Update handles structure
guidata(hObject, handles);

if handles.saveMov
    %disp('writeVideo is called..');
    currframe = getframe(hObject);%gcf is not appropriate??
    writeVideo(handles.mvobj, currframe);
end


% Utility function to show the next available frame
function ShowNextFrame(handleTimer, eventTimer)
hObject = handleTimer.UserData;
handles = guidata(hObject); 
iFrame = handles.CurrentFrame + 1;
if iFrame>handles.nFrames, 
    iFrame = 1; 
end
ShowFrame(hObject, handles, iFrame);

 

% Utility function to add a ROI
function AddOneROI(hObject,handles, roi)

lines = findall(roi,'type','line');%DS 2015/4/9
set(lines,'linewidth',1);%DS 2015/4/9
set(lines,'marker','.');%DS 2015/4/9

iRoi = length(handles.TimePlots)+1;

%DS 2015/4/9 new color order 
ColorList = ...  
    [0    0.4470    0.7410;...
    0.8500    0.3250    0.0980;...
    0.9290    0.6940    0.1250;...
    0.4940    0.1840    0.5560;...
    0.4660    0.6740    0.1880;...
    0.3010    0.7450    0.9330;...
    0.6350    0.0780    0.1840];
TheColor = ColorList(mod(iRoi-1,7)+1,:);

    
% ColorList = 'rcbmg';%NG ... hard to discriminate b and g
% TheColor = ColorList(mod(iRoi-1,5)+1);

handles.ROIs(iRoi) = roi;
handles.ROIs(iRoi).setColor(TheColor);

% add a default "delete" entry to the menu
hchild = get(roi, 'Children');
hhh = get(hchild(1),'UIContextMenu');
hhhh = get(hhh, 'Children');
set(hhhh(1),'visible','off');
%hchild(1).UIContextMenu.Children(1).delete %2015/4/10

hcmenu = get(hchild(1),'UIContextMenu');
itemnew = uimenu( hcmenu, ...
    'Label', 'Delete', ...
    'Callback', @(hObject,eventdata)MovieInspector('DelOneROI',hObject,guidata(hObject)),...
    'UserData', iRoi ); 

% I don't know how to deal with changes in color, so I am disabling it
foo = findobj(hcmenu,'Tag','set color cmenu item');
set(foo,'Enable','Off');

% we need to deal with:
% Changes in position
%   addNewPositionCallback(h,@(p) title(mat2str(p,3)));

handles.TimePlots(iRoi) = plot( handles.axTimeCourse, ...
        handles.tt, zeros(handles.nFrames,1),'Color', TheColor ); 

% % BUGGY clicking on it updates the time point:
% set(handles.TimePlots(iRoi),'ButtonDownFcn',...
%      @(hObject,eventdata)MovieInspector('axTimeCourse_ButtonDownFcn',hObject,eventdata,guidata(hObject)) );


addNewPositionCallback(roi,@(p) MovieInspector('ChangeOneROI',hObject,handles,iRoi,p));

% update the time course
ChangeOneROI(hObject,handles,iRoi,[])

% Update handles structure
guidata(hObject, handles);
 
% Utility function to reshape a ROI
function ChangeOneROI(hObject,handles,iRoi,p)

roi = handles.ROIs(iRoi);


% update the time course
coords = round(roi.getPosition);
coords(1) = max(coords(1),1);
coords(2) = max(coords(2),1);
coords(3) = min( coords(1)+coords(3), size(handles.mmm,2) ) - coords(1);
coords(4) = min( coords(2)+coords(4), size(handles.mmm,1) ) - coords(2);

% added by MC 2015-02-06 to ensure rectangle stays inside image
% (but oddly it does not seem to work)
% handles.ROIs(iRoi).setPosition(coords); 

AveragingRows = nanmean(  handles.mmm(coords(2)+[0:coords(4)],:,:),1);
AveragingCols = nanmean(AveragingRows(:,coords(1)+[0:coords(3)],:),2);
TimeCourse = squeeze(AveragingCols);

% ROImask = createMask(roi);
% [ii,jj,~]=find(ROImask);
% TimeCourse = squeeze(nanmean(nanmean(handles.mmm(ii,jj,:),1),2));
set(handles.TimePlots(iRoi),'YData',TimeCourse);

% Utility function to remove a ROI
function DelOneROI(hObject,handles)

iRoi = get(hObject,'UserData');

delete(handles.TimePlots(iRoi)); handles.TimePlots(iRoi) = [];
delete(handles.ROIs     (iRoi)); handles.ROIs     (iRoi) = [];

% Update handles structure
guidata(handles.figMovieInspector, handles);


% --- Executes on mouse press over axes background.
function axTimeCourse_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axTimeCourse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axTimeCourse = hObject;  % MC
CurrentPoint = get(axTimeCourse,'CurrentPoint'); % oddly, 2 by 3
DesiredTime = CurrentPoint(1);
if DesiredTime < handles.tt(1), %DS modified on 3/5/15
    % for now this does nothing to change the actual ylims!!
    fprintf('Calculating histograms...');
    ylims = double(get(handles.axTimeCourse,'ylim'));
    figure; clf; ax = gridplot(2,1);
    axes(ax(1));
    ff = isfinite(handles.mmm(:));
    hist(handles.mmm(ff)); 
    text( ylims(1), 0, '\downarrow', 'hori', 'center', 'vert', 'bottom', 'fontsize', 20, 'color', 'r' )
    text( ylims(2), 0, '\downarrow', 'hori', 'center', 'vert', 'bottom', 'fontsize', 20, 'color', 'g' )
    axes(ax(2));
    xx = ylims+[-1 1]*diff(ylims)/10;
    hist(handles.mmm(ff),linspace(xx(1),xx(2),20));
    text( ylims(1), 0, '\downarrow', 'hori', 'center', 'vert', 'bottom', 'fontsize', 20, 'color', 'r' )
    text( ylims(2), 0, '\downarrow', 'hori', 'center', 'vert', 'bottom', 'fontsize', 20, 'color', 'g' )
    set(ax(2),'xlim',xx);
    fprintf('done.\n');

    return; 
end
[~, iFrame] = min(abs(handles.tt - DesiredTime));
ShowFrame(hObject, handles,iFrame);

% --- Executes on mouse press over axes background.
function axMovie_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axMovie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

CurrentPoint = get(handles.axMovie,'CurrentPoint'); % oddly, 2 by 3
Coords = round(CurrentPoint(1,[2 1]));
fprintf('Selected ROI around row %d, column %d\n',Coords);
Sz = 2*ceil( min(size(handles.mmm,1),size(handles.mmm,2))/4 ); % round and divisible by 2
roi = imrect (handles.axMovie, [Coords(2)-Sz/2,Coords(1)-Sz/2,Sz,Sz]);

AddOneROI(hObject,handles,roi);

    
% --- Executes on mouse press over figure background.
function figMovieInspector_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to figMovieInspector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fprintf('Click on the image to select a ROI.\nClick on the timeline to select a time point.\n');


% --- Executes when user attempts to close figMovieInspector.
function figMovieInspector_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figMovieInspector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

delete(handles.Timer);

if handles.saveMov
    close(handles.mvobj);
end

% Hint: delete(hObject) closes the figure
delete(hObject);


% --- Executes on button press in pbSave.
function pbSave_Callback(hObject, eventdata, handles)
% hObject    handle to pbSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~handles.saveMov
    
    handles = guidata(hObject); 
    handles.saveMov = 1;
    guidata(hObject, handles);
    
disp('Recording mode ON. Press ">" to start, "||" to stop ');

%     if handles.Timer.Running, stop(handles.Timer); end
%     set(handles.Timer,'Period',0.1);
%     handles.CurrentFrame = 1;
    
    open(handles.mvobj);
    
    %start(handles.Timer);

end


% --- Executes on slider movement.
function sliderSpeed_Callback(hObject, eventdata, handles)
% hObject    handle to sliderSpeed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
timerPeriod = round(1000* median(diff(handles.tt))/10^(get(hObject, 'Value')))/1000;
set(handles.text1, 'String', ['x ' num2str(median(diff(handles.tt))/timerPeriod,2)]); 

if strcmp(handles.Timer.Running,'off')
    set(handles.Timer,'Period', timerPeriod);
elseif handles.Timer.Running, stop(handles.Timer);
    set(handles.Timer,'Period', timerPeriod);
    set(handles.pbPlayPause, 'String', '||');
    %handles.figMovieInspector.Children(length(handles.figMovieInspector.Children) - 2).String = '||';
    start(handles.Timer);
end
% % --- Executes during object creation, after setting all properties.
function sliderSpeed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderSpeed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% set relative speed of movie, compared to recorded movie

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

set(hObject, 'Min', log10(0.1));
set(hObject, 'Max', log10(50));
set(hObject, 'Value', log10(1));
