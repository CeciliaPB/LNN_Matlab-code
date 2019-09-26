function MClustCutterRedrawAxes(varargin)

% Fork
if ~ischar(varargin{1})
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
        'Position',MClust_CHDrawingAxisWindow_Pos);   % I do not have 'MClustCutterKeyPress'
else
    % figure already exists -- select it
    figure(drawingFigHandle);
end

% have to a complete redraw
if ~full
    curAxis = axis;
end
cla;
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
                if usefastplot 
                    h = fastplot(MClust_CurrentFeatureData(f,1), MClust_CurrentFeatureData(f,2), marker);
                else
                    h = plot(MClust_CurrentFeatureData(f,1), MClust_CurrentFeatureData(f,2), marker);
                end
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
            if usefastplot
                h = fastplot(MClust_CurrentFeatureData(f,1), MClust_CurrentFeatureData(f,2), marker);
            else
                h = plot(MClust_CurrentFeatureData(f,1), MClust_CurrentFeatureData(f,2), marker);
            end
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
if ~usefastplot
    zoom on
end

contourWindow = findobj('Type', 'figure', 'Tag', 'ContourWindow');
if ~isempty(contourWindow)
    mkContours(drawingFigHandle, 'figHandle', contourWindow);
end
figure(drawingFigHandle);

% Store variables
setappdata(gcf,'plot_input_arguments',varargin)
setappdata(gcf,'control_window_handle',figHandle)
setappdata(gcf,'xrange',[min(MClust_CurrentFeatureData(:,1)) max(MClust_CurrentFeatureData(:, 1))+0.0001])
setappdata(gcf,'yrange',[min(MClust_CurrentFeatureData(:,2)) max(MClust_CurrentFeatureData(:, 2))+0.0001])

% -------------------------------------------------------------------------
function zoom_main(figHandle, varargin)

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

cla;
hold on;
for iC = 0:nClust
    if ~MClust_Hide(iC+1)
        HideClusterHandle = findobj(figHandle, 'UserData', iC, 'Tag', 'HideCluster');
        if iC == 0
            if MClust_UnaccountedForOnly
                MClust_ClusterIndex = ProcessClusters(MClust_CurrentFeatureData, MClust_Clusters);
                f = (MClust_ClusterIndex == 0);
                figure(drawingFigHandle);
                if usefastplot 
                    h = fastplot(MClust_CurrentFeatureData(f,1), MClust_CurrentFeatureData(f,2), marker);
                else
                    h = plot(MClust_CurrentFeatureData(f,1), MClust_CurrentFeatureData(f,2), marker);
                end
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
            if usefastplot
                h = fastplot(MClust_CurrentFeatureData(f,1), MClust_CurrentFeatureData(f,2), marker);
            else
                h = plot(MClust_CurrentFeatureData(f,1), MClust_CurrentFeatureData(f,2), marker);
            end
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
xlabel(MClust_CurrentFeatureNames{1},'interpreter','none');
ylabel(MClust_CurrentFeatureNames{2},'interpreter','none');

contourWindow = findobj('Type', 'figure', 'Tag', 'ContourWindow');
if ~isempty(contourWindow)
    mkContours(drawingFigHandle, 'figHandle', contourWindow);
end

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
% xpixels = pos(3) * 2;     % no rendering 'stripes', but slower
% ypixels = pos(4) * 2;

mnx = min(x0);
mxx = max(x0);
mny = min(y0);
mxy = max(y0);
x2 = round((x0-mnx)/(mxx-mnx)*(xpixels-1)) + 1;
y2 = round((y0-mny)/(mxy-mny)*(ypixels-1)) + 1;
u = unique(x2*100000+y2);
y3 = mod(u,100000);
x3 = (u - y3) / 100000;
x4 = (x3 / xpixels) * (mxx - mnx) + mnx;
y4 = (y3 / ypixels) * (mxy - mny) + mny;

h = plot(x4,y4,marker);

% Restore axis units property
set(gca,'Unit',old_units)
% set(gcf,'Renderer','openGL')    % openGL renders somewhat nicer than zbuffers

% Set ResizeFcn
rsf = 'MClustCutterRedrawAxes(''figure_ResizeFcn'',gcf)';
set(gcf,'ResizeFcn',rsf)

% Set ButtonDownFcn
bdf = 'MClustCutterRedrawAxes(''figure_ButtonDownFcn'',gcf)';
set(gca,'ButtonDownFcn',bdf)
set(allchild(gca),'ButtonDownFcn',bdf)

% -------------------------------------------------------------------------
function figure_ResizeFcn(hObj)     %#ok<DEFNU>

% Get input arguments used to create the figure
plargs = getappdata(gcf,'plot_input_arguments');
figHandle = getappdata(gcf,'control_window_handle');

% Redraw
main(figHandle,plargs{:});

% -------------------------------------------------------------------------
function figure_ButtonDownFcn(hObj)     %#ok<DEFNU>

% Get variables
A = gca;
xx = get(A,'XLim');
yy = get(A,'YLim');
xrange = getappdata(gcf,'xrange');
yrange = getappdata(gcf,'yrange');

% Set axis
seltyp = get(hObj,'SelectionType');
switch seltyp
case 'normal'   % zoom in
    point1 = get(A,'CurrentPoint'); % button down detected
    rbbox;
    point2 = get(A,'CurrentPoint');    % button up detected
    point1 = point1(1,1:2);              % extract x and y
    point2 = point2(1,1:2);
    if isequal(point1,point2)
        xx2 = abs(xx(2)-xx(1)) / 4;
        yy2 = abs(yy(2)-yy(1)) / 4;
        xx3(1) = point1(1) - xx2;
        xx3(2) = point1(1) + xx2;
        yy3(1) = point1(2) - yy2;
        yy3(2) = point1(2) + yy2;
        xnew(1) = max(xrange(1),xx3(1));
        ynew(1) = max(yrange(1),yy3(1));
        xnew(2) = min(xrange(2),xx3(2));
        ynew(2) = min(yrange(2),yy3(2));
    else
        xnew = [min([point1(1) point2(1)]) max([point1(1) point2(1)])];
        ynew = [min([point1(2) point2(2)]) max([point1(2) point2(2)])];
    end

case 'open'   % set default
    xnew = xrange;
    ynew = yrange;
    
case 'extend'   % zoom out
    point = get(A,'CurrentPoint');   % button down detected
    point = point(1,1:2);
    xx2 = abs(xx(2)-xx(1));
    yy2 = abs(yy(2)-yy(1));
    xx3(1) = point(1) - xx2;
    xx3(2) = point(1) + xx2;
    yy3(1) = point(2) - yy2;
    yy3(2) = point(2) + yy2;
    xnew(1) = max(xrange(1),xx3(1));
    ynew(1) = max(yrange(1),yy3(1));
    xnew(2) = min(xrange(2),xx3(2));
    ynew(2) = min(yrange(2),yy3(2));
        
otherwise
    return
end

if xnew(1) >= xnew(2)    % problems can occur when too zoomed in to resolve
    xnew(2) = xnew(1) + eps;
end
if ynew(1) >= ynew(2)
    ynew(2) = ynew(1) + eps;
end
set(A,'XLim',xnew);
set(A,'YLim',ynew);

% Get input arguments used to create the figure
plargs = getappdata(gcf,'plot_input_arguments');
figHandle = getappdata(gcf,'control_window_handle');

% Redraw
zoom_main(figHandle,plargs{:},'full',0);

% global MClust_CurrentFeatureData
% h = fastplot(MClust_CurrentFeatureData(:,1),MClust_CurrentFeatureData(:,2),'.')