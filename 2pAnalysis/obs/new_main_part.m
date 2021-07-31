function pushbutton2_Callback(hObject, eventdata, h)
% variance explained mask
h.dat.cl.vmap = 'unit';
set_maskBcolor(h, 1)

redraw_figure(h);
guidata(hObject,h);

% function Q33_Callback(hObject, eventdata, h)
% iy = 3; ix = 3;
% quadrant(hObject, h, iy, ix);
% paint_quadbutton(h, iy, ix);

% function paint_quadbutton(h, iy, ix)
% set(h.full, 'BackgroundColor', .92 * [1 1 1])
% 
% for j = 1:3
%     for i = 1:3
%         if h.quadvalue(j,i)==1
%             set(h.(sprintf('Q%d%d', j,i)), 'BackgroundColor', [.4 .4 .4]);
%         end
%     end
% end
% set(h.(sprintf('Q%d%d', iy,ix)), 'BackgroundColor', [1 0 0]);

% need this??
% function full_Callback(hObject, eventdata, h)
% h.dat.ylim = [0 h.dat.cl.Ly];
% h.dat.xlim = [0 h.dat.cl.Lx];
% if h.dat.map==1
%     redraw_figure(h);
% else
%     redraw_meanimg(h);
% end
% set(h.full, 'BackgroundColor', [1 0 0]);
% for i = 1:3
%     for j = 1:3
%         if h.(sprintf('Q%d%d', j,i)).BackgroundColor(1) >.99
%             set(h.(sprintf('Q%d%d', j,i)), 'BackgroundColor', [.4 .4 .4])
%         end
%     end
% end
% guidata(hObject,h);

% function quadrant(hObject, h, iy, ix)
% h.dat.ylim = [h.dat.figure.y0all(iy) h.dat.figure.y1all(iy+1)];
% h.dat.xlim = [h.dat.figure.x0all(ix) h.dat.figure.x1all(ix+1)];
% h.quadvalue(iy, ix) = 1;
% 
% guidata(hObject,h);
% 
% if h.dat.map==1
%     redraw_figure(h);
% else
%     redraw_meanimg(h);
% end

%% short cut key
function figure1_WindowKeyPressFcn(hObject, eventdata, h)
switch eventdata.Key  
    case 'q' %roi random color
        pushbutton87_Callback(hObject, eventdata, h);
    case 'e' %mean image
        pushbutton89_Callback(hObject, eventdata, h);
%     case 'z'
%         pushbutton102_Callback(hObject, eventdata, h);
end


%% CELL CLICKING!! 
function figure1_WindowButtonDownFcn(hObject, eventdata, h)
% house the clicked cell id as h.dat.Fchosen
% other changed values in h
% h.text54 information on the clicked cell is displayed at bottom left
% other variables used h.dat
% res.iclust
% cl.Lx
% cl.Ly

z = round(eventdata.Source.CurrentAxes.CurrentPoint(1,:));
x = round(z(1));
y  = round(z(2));

if x>=1 && y>=1 && x<=h.dat.cl.Lx && y<=h.dat.cl.Ly && h.dat.res.iclust(y,x)>0
    h.dat.ichosen = h.dat.res.iclust(y, x);
    ichosen = h.dat.ichosen;
     
    redraw_figure_test(h);
    guidata(hObject,h);
    
    str = sprintf('cell #%d\n',ichosen);
    
    %will be deleted
% %     labels = [h.statLabels(2:end), {'iscell'}, {'redcell'},{'redprob'}];
% %     for j =1:length(labels)
% %         if isfield(h.dat.stat, labels{j})
% %             sl = eval(sprintf('h.dat.stat(ichosen).%s', labels{j}));
% %             strnew = sprintf('%s = %2.2f \n', labels{j}, sl);
% %             str = cat(2, str, strnew);
% %         end
% %     end
  
%will be recovered
 %   set(h.text54,'String', str);
end




%------------------------ BACKGROUND --------------------------------%

% function set_Bcolor(h, ih)
% pb = [87 103 89 90 92];
%
% for j = 1:length(pb)
%     if j==ih
%         set(h.(sprintf('pushbutton%d', pb(ih))), 'BackgroundColor', [1 0 0]);
%     else
%         if h.(sprintf('pushbutton%d', pb(j))).BackgroundColor(1)>.99
%             set(h.(sprintf('pushbutton%d', pb(j))), 'BackgroundColor', .94 * [1 1 1]);
%         end
%     end
% end

% roi random color
function pushbutton87_Callback(hObject, eventdata, h)
% h.dat.map = 1;
redraw_figure(h);
%set_Bcolor(h, 1);
guidata(hObject,h);

% mean img
function pushbutton89_Callback(hObject, eventdata, h)
%h.dat.map = 2;
redraw_meanimg_test(h);
%set_Bcolor(h, 3);
guidata(hObject,h);

% % --- meanimg
% function pushbutton102_Callback(hObject, eventdata, h)
% h=meanimgFig(h);
% guidata(hObject,h);
% 
% function h=meanimgFig(h)
% hval0            = [h.dat.stat.mimgProjAbs];
% hval             = hval0;
% hval             = hval - min(hval);
% hval             = hval / nanmean(hval);
% hval             = .6 * (1 - min(2*nanstd(hval)+1, hval)/(nanstd(hval)*2 + 1)) + .15;
% h.dat.cl.rands   = hval;
% h.dat.cl.cmap    = [hval0(:) hval(:)];
% I = redraw_figure(h);
% h.dat.cl.meanimgFig = I;
% set_maskCcolor(h, 4);

%------------------------ BACKGROUND --------------------------------%

% function set_Bcolor(h, ih)
% pb = [87 103 89 90 92];
%
% for j = 1:length(pb)
%     if j==ih
%         set(h.(sprintf('pushbutton%d', pb(ih))), 'BackgroundColor', [1 0 0]);
%     else
%         if h.(sprintf('pushbutton%d', pb(j))).BackgroundColor(1)>.99
%             set(h.(sprintf('pushbutton%d', pb(j))), 'BackgroundColor', .94 * [1 1 1]);
%         end
%     end
% end