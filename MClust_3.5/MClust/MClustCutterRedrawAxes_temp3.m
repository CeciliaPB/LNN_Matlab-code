function MClustCutterRedrawAxes(varargin)

% Fork
if isa(varargin{1},'figure')    % call the zoom
    axes_zoom(varargin{:})
elseif ~ischar(varargin{1})
    main(varargin{1},varargin{2:end})
else    % invoke named subfunction or callback
	feval(varargin{:});
end

% -------------------------------------------------------------------------
function main(figHandle, varargin)

usefastplot = 1;

global MClust_Clusters MClust_Colors MClust_Hide MClust_UnaccountedForOnly 
global MClust_ClusterIndex MClust_FeatureData 

global MClust_CurrentFeatures % used to keep track of which features are currently in memory
global MClust_CurrentFeatureNames % 
global MClust_CurrentFeatureData
global MClust_FeatureNames % names of features
global MClust_FeatureSources % <filenames, number pairs> for finding features in fd files

global MClust_ClusterCutWindow_Marker
global MClust_ClusterCutWindow_MarkerSize
global MClust_CHDrawingAxisWindow_Pos;

% -- get variables
full = 0;
extract_varargin;

nClust = length(MClust_Clusters);

drawingFigHandle = findobj('Type', 'figure', 'Tag', 'CHDrawingAxisWindow');  % figure to draw in

xdimHandle = findobj(figHandle, 'Tag', 'xdim');
xdim = get(xdimHandle, 'Value');           % x dimemsion to plot
ydimHandle = findobj(figHandle, 'Tag', 'ydim');  
ydim = get(ydimHandle, 'Value');           % y dimension to plot
markerHandle = findobj(figHandle, 'Tag', 'PlotMarker');
markerString = get(markerHandle, 'String');
markerValue = get(markerHandle, 'Value');
MClust_ClusterCutWindow_Marker = markerValue;
marker = markerString{markerValue};
markerSizeHandle = findobj(figHandle, 'Tag', 'PlotMarkerSize');
markerSizeString = get(markerSizeHandle, 'String');
markerSizeValue = get(markerSizeHandle, 'Value');
MClust_ClusterCutWindow_MarkerSize = markerSizeValue;
markerSize = str2double(markerSizeString{markerSizeValue});

% converted back to work by disk access (ADR 2008 - turns out this is faster
% get xdim
if (MClust_CurrentFeatures(1) ~= xdim)
    temp = load(MClust_FeatureSources{xdim,1}, '-mat', 'FeatureData');
    MClust_CurrentFeatureData(:,1) = temp.FeatureData(:,MClust_FeatureSources{xdim,2});
    MClust_CurrentFeatures(1) = xdim;
    MClust_CurrentFeatureNames{1} = MClust_FeatureNames{xdim};
end
% get ydim
if (MClust_CurrentFeatures(2) ~= ydim)
    temp = load(MClust_FeatureSources{ydim,1}, '-mat', 'FeatureData');
    MClust_CurrentFeatureData(:,2) = temp.FeatureData(:,MClust_FeatureSources{ydim,2});
    MClust_CurrentFeatures(2) = ydim;
    MClust_CurrentFeatureNames{2} = MClust_FeatureNames{ydim};
end

if isempty(drawingFigHandle)
    % create new drawing figure
    drawingFigHandle = ...
        figure('Name', 'Cluster Cutting Window',...
        'NumberTitle', 'off', ...
        'Tag', 'CHDrawingAxisWindow', ...
        'KeyPressFcn', 'MClustCutterKeyPress','Position',MClust_CHDrawingAxisWindow_Pos);
else
    % figure already exists -- select it
    figure(drawingFigHandle);
end

% have to a complete redraw
if ~full
    curAxis = axis;
end
clf;
hold on;
if full %%% Added by JCJ Aug 2007 to stabilize redraw of display
    set(gca, 'XLim', [min(MClust_CurrentFeatureData(:,1)) max(MClust_CurrentFeatureData(:,1))+0.0001]);
    set(gca, 'YLim', [min(MClust_CurrentFeatureData(:,2)) max(MClust_CurrentFeatureData(:,2))+0.0001]);	
else
    axis(curAxis);
end
for iC = 0:nClust
    if ~MClust_Hide(iC+1)
        HideClusterHandle = findobj(figHandle, 'UserData', iC, 'Tag', 'HideCluster');
        if iC == 0
            if MClust_UnaccountedForOnly
                MClust_ClusterIndex = ProcessClusters(MClust_CurrentFeatureData, MClust_Clusters);
                f = (MClust_ClusterIndex == 0);
                 figure(drawingFigHandle);
                h = plot(MClust_CurrentFeatureData(f,1), MClust_CurrentFeatureData(f,2), marker);
            else
                figure(drawingFigHandle);
                allf = 1:size(MClust_CurrentFeatureData,1);
                sbst = allf(1:end);
                if usefastplot
                    h = fastplot(MClust_CurrentFeatureData(sbst,1), MClust_CurrentFeatureData(sbst,2), marker);
                else
                    h = plot(MClust_CurrentFeatureData(sbst,1), MClust_CurrentFeatureData(sbst,2), marker);
                end
            end
        else         
            [f,MClust_Clusters{iC}] = FindInCluster(MClust_Clusters{iC});
            if isempty(f) && ~isempty(HideClusterHandle)
                set(HideClusterHandle, 'Enable', 'off');
            else 
                set(HideClusterHandle, 'Enable', 'on');
            end
            figure(drawingFigHandle);
            h = plot(MClust_CurrentFeatureData(f,1), MClust_CurrentFeatureData(f,2), marker);
        end
        set(h, 'Color', MClust_Colors(iC+1,:));
        set(h, 'Tag', 'ClusterLine', 'UserData', iC);
        set(h, 'MarkerSize', markerSize);
		if iC > 0
			try
				DrawOnAxis(MClust_Clusters{iC}, xdim, ydim, MClust_Colors(iC+1,:), gca); 
			end
		end
    end
end
figure(drawingFigHandle);
if full
    set(gca, 'XLim', [min(MClust_CurrentFeatureData(:,1)) max(MClust_CurrentFeatureData(:, 1))+0.0001]);
    set(gca, 'YLim', [min(MClust_CurrentFeatureData(:,2)) max(MClust_CurrentFeatureData(:, 2))+0.0001]);
else
    axis(curAxis);
end
xlabel(MClust_CurrentFeatureNames{1},'interpreter','none');
ylabel(MClust_CurrentFeatureNames{2},'interpreter','none');
% zoom on

contourWindow = findobj('Type', 'figure', 'Tag', 'ContourWindow');
if ~isempty(contourWindow)
    mkContours(drawingFigHandle, 'figHandle', contourWindow);
end
figure(drawingFigHandle);

% Store variables
setappdata(gcf,'plot_input_arguments',varargin)
setappdata(gcf,'control_window_handle',figHandle)


% -------------------------------------------------------------------------
function h = fastplot(x,y,marker)

% Restrict according to the axis limits
xl = xlim;
yl = ylim;
rinx = x > xl(1) & x < xl(2) & y > yl(1) & y < yl(2);
x0 = x(rinx);
y0 = y(rinx);

% Plot
old_units = get(gca,'Units');
set(gca,'Units','pixels')
pos = get(gca,'Position');
xpixels = pos(3);
ypixels = pos(4);

mnx = min(x0);
mxx = max(x0);
mny = min(y0);
mxy = max(y0);
x2 = round((x0-mnx)/(mxx-mnx)*xpixels);
y2 = round((y0-mny)/(mxy-mny)*ypixels);
u = unique(x2*100000+y2);
y3 = mod(u,100000);
x3 = (u - y3) / 100000;
x4 = (x3 / xpixels) * (mxx - mnx) + mnx;
y4 = (y3 / ypixels) * (mxy - mny) + mny;

h = plot(x4,y4,marker);

% Restore axis units property
set(gca,'Unit',old_units)

% Set ResizeFcn
rsf = 'MClustCutterRedrawAxes(''figure_ResizeFcn'',gcf)';
set(gcf,'ResizeFcn',rsf)

% Zoom
z = zoom;
zcb = 'MClustCutterRedrawAxes(''axes_zoom'',gcf)';
set(z,'ActionPostCallback',@MClustCutterRedrawAxes);
set(z,'Enable','on')

% Set ButtonDownFcn
% bdf = 'MClustCutterRedrawAxes(''figure_ButtonDownFcn'',gcf)';
% set(gca,'ButtonDownFcn',bdf)
% set(allchild(gca),'ButtonDownFcn',bdf)

% -------------------------------------------------------------------------
function figure_ResizeFcn(hObj)     %#ok<DEFNU>

% Get input arguments used to create the figure
plargs = getappdata(gcf,'plot_input_arguments');
figHandle = getappdata(gcf,'control_window_handle');

% Redraw
main(figHandle,plargs{:});

% -------------------------------------------------------------------------
function axes_zoom(hObj,ev)

% Get input arguments used to create the figure
plargs = getappdata(gcf,'plot_input_arguments');
figHandle = getappdata(gcf,'control_window_handle');

% Redraw
main(figHandle,plargs{:},'full',0);

% -------------------------------------------------------------------------
function figure_ButtonDownFcn(hObj)     %#ok<DEFNU>

% Get variables
A = gca;
binraster = getappdata(hObj,'binraster');
time = getappdata(hObj,'time');
tno = getappdata(hObj,'tno');
tno0 = getappdata(hObj,'tno0');
time0 = getappdata(hObj,'time0');
rtime0 = getappdata(hObj,'realtime0');

% Set axis
seltyp = get(hObj,'SelectionType');
switch seltyp
case 'normal'   % zoom in
    point1 = get(A,'CurrentPoint'); % button down detected
    rbbox;
    point2 = get(A,'CurrentPoint');    % button up detected
    point1 = point1(1,1:2);              % extract x and y
    point2 = point2(1,1:2);
    point1(1) = linterp(rtime0,time0,point1(1));
    point2(1) = linterp(rtime0,time0,point2(1));
    if isequal(point1,point2)
        xx = [time(1) time(end)];
        yy = [tno(1) tno(end)];
        xx2 = (abs(xx(2) - xx(1))) / 4;
        yy2 = (abs(yy(2) - yy(1))) / 4;
        xx3(1) = point1(1) - xx2;
        xx3(2) = point1(1) + xx2;
        yy3(1) = point1(2) - yy2;
        yy3(2) = point1(2) + yy2;
        if time0(1) < xx3(1) && time0(end) > xx3(2)
            xlim_data = xx3;
        elseif time0(1) > xx3(1)
            xx_new(1) = time(1);
            xx_new(2) = time(1) + (2 * xx2);
            xlim_data = xx_new;
        elseif time0(end) < xx3(2)
            xx_new(1) = time0(end) - (2 * xx2);
            xx_new(2) = time0(end);
            xlim_data = xx_new;
        end
        if tno0(1) < yy3(1) && tno0(end) > yy3(2)
            ylim_data = yy3;
        elseif tno0(1) > yy3(1)
            yy_new(1) = tno0(1);
            yy_new(2) = tno0(1) + (2 * yy2);
            ylim_data = yy_new;
        elseif tno0(end) < yy3(2)
            yy_new(1) = tno0(end) - (2 * yy2);
            yy_new(2) = tno0(end);
            ylim_data = yy_new;
        end
    else
        xlim_data = [min([point1(1) point2(1)]) max([point1(1) point2(1)])];
        xlim_data(1) = max(xlim_data(1),time0(1));
        xlim_data(2) = min(xlim_data(2),time0(end));
        ylim_data = [min([point1(2) point2(2)]) max([point1(2) point2(2)])];
        ylim_data(1) = max(ylim_data(1),tno0(1));
        ylim_data(2) = min(ylim_data(2),tno0(end));
    end

case 'open'   % set default
    xlim_data = [time0(1) time0(end)];
    ylim_data = [tno0(1) tno0(end)];
    
case 'extend'   % zoom out
    point = get(A,'CurrentPoint'); % button down detected
    point = point(1,1:2);
    point(1) = linterp(rtime0,time0,point(1));
    xx = [time(1) time(end)];
    yy = [tno(1) tno(end)];
    xx2 = abs(xx(2) - xx(1));
    yy2 = abs(yy(2) - yy(1));
    if xx2 > (time0(end) - time0(1)) / 2,
        xx2 = (time0(end) - time0(1)) / 2;
    end
    if yy2 > (abs(tno0(end) - tno0(1))) / 2,
        yy2 = (abs(tno0(end) - tno0(1))) / 2;
    end
    xx3(1) = point(1) - xx2;
    xx3(2) = point(1) + xx2;
    yy3(1) = point(2) - yy2;
    yy3(2) = point(2) + yy2;
    if time0(1) < xx3(1) && time0(end) > xx3(2)
        xlim_data = xx3;
    elseif time0(1) > xx3(1)
        xx_new(1) = time0(1);
        xx_new(2) = time0(1) + (2 * xx2);
        xlim_data = xx_new;
    elseif time0(end) < xx3(2)
        xx_new(1) = time0(end) - (2 * xx2);
        xx_new(2) = time0(end);
        xlim_data = xx_new;
    end
    if tno0(1) < yy3(1) && tno0(end) > yy3(2)
        ylim_data = yy3;
    elseif tno0(1) > yy3(1)
        yy_new(1) = tno0(1);
        yy_new(2) = tno0(1) + (2 * yy2);
        ylim_data = yy_new;
    elseif tno0(end) < yy3(2)
        yy_new(1) = tno0(end) - (2 * yy2);
        yy_new(2) = tno0(end);
        ylim_data = yy_new;
    end
otherwise
    return
end

% Redraw
time2 = round(xlim_data(1)):round(xlim_data(2));
tno2 = round(ylim_data(1)):round(ylim_data(2));
rtime2 = rtime0(time2(1):time2(end));
binraster2 = binraster(tno2,time2);
setappdata(hObj,'realtime',rtime2);
setappdata(hObj,'time',time2);
setappdata(hObj,'tno',tno2);
main(hObj,binraster2);

% -------------------------------------------------------------------------
function yi = linterp(x,y,xi)
%LINTERP   1D linear interpolation.
%   YI = LINTERP(X,Y,XI) interpolates Y = f(XI), where f is given in (X,Y).
%
%   See also INTERP1Q.

% Interpolation
lx = length(xi);
lex = length(x);
yi = zeros(1,lx);
for k = 1:lx    % linear interpolation
    inx1 = find(x<=xi(k),1,'last');
    inx1(inx1==lex) = inx1(inx1==lex) - 1;
    inx2 = inx1 + 1;
    yi(k) = y(inx1) + (y(inx2) - y(inx1)) * (xi(k) - x(inx1)) / (x(inx2) - x(inx1));
end