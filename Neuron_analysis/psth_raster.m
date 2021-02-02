function [rastmat, timevec] = psth_raster(trialspx,pre,post,varargin)

% Generates a matrix for raster plot from cell array 'trialspx' generated
% by main function.
%
% INPUTS 
% - trialspx: vector with timestamps of spike events (units as recorded)
% - pre: time before trigger to include in psth. Default 1s. 
% - post: time after trigger to include in psth. Default 1s.  
% Varargin 
% - 'chart': generate a plot. Default = 0, no plot will be generated. If 1,
% a plot will be generated
%
% OUTPUT
% - rastmat: peri-stimulus time histogram
% - timevec: peri-stimulus time histogram up to the first after the trigger
% 
% Examples 
% [rastmat,timevec] = mraster(trialspx,1000,1000,'chart')
%
% -------------------------------------------------------------------------
% Cecília Pardo-Bellver, 2019
% Laboratory of Network Neurophysiology
% Institute of Experimantal Medicine, Hungary.
%
% Based on mraster by Maik C. Stttgen, Feb 2011
% -------------------------------------------------------------------------

% Read varargin
chart=0;
if nargin>3
  if strcmp(varargin{1},'chart')
    chart = 1;
  end
end

% preallocate
rastmat = zeros(numel(trialspx),pre+1+post);
timevec = -pre:1:post;  

% generate raster
for ii = 1:numel(trialspx)
  rastmat(ii,trialspx{ii}+pre+1) = 1;
end

% plot raster
if chart == 1
  figure('name','Peri-stimulus time histogram','units','normalized','position',[0.3 0.4 0.4 0.2])
  for ii = 1:numel(trialspx)
    plot(timevec(rastmat(ii,:)~=0),rastmat(ii,rastmat(ii,:)~=0)*ii,'k.','MarkerSize',4);
    hold on; 
  end
  axis([-pre+10 post+10 0.5 numel(trialspx)+0.5])
  xlabel('Time (ms)'),ylabel('trials')
end

