function MCC = createClusterObject(Name,TimeStamps)
%CREATECLUSTEROBJECT   Create a cluster object from time stamps.
%   MCC = CREATECLUSTEROBJECT(NAME,TIMESTAMPS) creates an MClust cluster
%   object (class name: 'mccluster') with the given name and time stamp.
%
%   See also MCCCLUSTER.

% Convert to column vector
TimeStamps = TimeStamps(:);

% Define properties
MCC.name = Name;
MCC.xdimNames = {}; 
MCC.ydimNames = {}; 
MCC.xdimSources = {};
MCC.ydimSources = {};
MCC.cx = {}; 
MCC.cy = {};
MCC.AddFlag = [];
MCC.recalc = 0;
MCC.myPoints = TimeStamps;
MCC.myOrigPoints = TimeStamps;
MCC.ForbiddenPoints = [];

% Make it an object of the mccluster class
MCC = class(MCC, 'mccluster');