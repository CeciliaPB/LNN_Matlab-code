function [x,y] = DrawConvexHull

% [x,y] = DrawConvexHull
% allows the user to draw a convex hull on the current axis
% and returns the x,y points on that hull
%
% ADR 1998
% Status PROMOTED
% version V4.1
%
% RELEASED as part of MClust 2.0
% See standard disclaimer in Contents.m
% 
% ADR fixed to handle empty inputs

usemycode = 1;

if usemycode
    % Draw polygon
seltyp = get(handles.figure1,'SelectionType');  % click type
if isequal(seltyp,'normal')
    point1 = get(handles.axes1,'CurrentPoint'); % button down detected
    point1x = point1(1,1);
    point1y = point1(1,2);
    point0x = point1(1,1);
    point0y = point1(1,2);
    L = [];
    bp = 1;
    while bp
        bp = waitforbuttonpress;    % wait for mouse click
    end
    seltyp2 = get(handles.figure1,'SelectionType');  % click type
    while isequal(seltyp2,'normal')
        point2 = get(handles.axes1,'CurrentPoint'); % button down detected
        point2x = point2(1,1);
        point2y = point2(1,2);
        L(end+1) = line([point1x point2x],[point1y point2y],'Color','c','LineWidth',2);
        point0x = [point0x point2x];
        point0y = [point0y point2y];
        point1x = point2x;
        point1y = point2y;
        bp = 1;
        while bp
            bp = waitforbuttonpress;    % wait for mouse click
        end
        seltyp2 = get(handles.figure1,'SelectionType');  % click type
    end
    L(end+1) = line([point2x point0x(1)],[point2y point0y(1)],'Color','c','LineWidth',2);
end
x = point0x;
y = point0y;
% delete(L)   % delete lines
end

else
[x,y] = ginput;
if isempty(x) || isempty(y)
    return
end
end
k = convexhull(x,y);
x = x(k);
y = y(k);