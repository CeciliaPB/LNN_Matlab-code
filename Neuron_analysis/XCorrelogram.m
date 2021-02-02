function [Y,X] = XCorrelogram(n1, n2, varargin)

% This function cauculates the cross correlation between two neurons.
%
% INPUT: 
%   - n1: .mat file with thw timestamps of neuron 1. ex. 'GR1_1.mat'
%   - n2: .mat file with thw timestamps of neuron 2. ex. 'GR1_2.mat'
%   Varargin
%   - 'bins': bin size of choice in msec. Default 1.
%   - 'width': msec on each side of the xcorr plot. Default 250.
%   - 's': saves the figure in the specified file format (jpg, png,
%   svg, etc.). See saveas help for more info. Default NOT save.
%
% OUTPUT:
%   - Y: y dim of the plot
%   - X: x dim of the plot, in msec.
%   - without output arguments, the function plots the xcorr
% 
% Examples: 
% [Y, X] = XCorrelogram(n1, n2, 'bin', 1, 'width', 250, 's', 'jpg');
% XCorrelogram(n1, n2, 's', 'jpg');
% 
% -------------------------------------------------------------------------
% Cecília Pardo-Bellver, 2019
% Laboratory of Network Neurophysiology
% Institute of Experimantal Medicine, Hungary.
%
% Uses code from MClust, MEX-file CrossCorr.
% -------------------------------------------------------------------------

load(n1,'TS');
TS1 = TS;
load(n2,'TS');
TS2 = TS;

% Default params
bin = 1;
width = 250;
s = 0;
    
% Params, read from varargin
if nargin
  for i=1:size(varargin,2)
    switch varargin{i}
      case 'bins'
        bin = varargin{i+1};
      case 'width'
        width = varargin{i+1};   
      case 's'
        s = 1;
        type = varargin{i+1};
      otherwise  
    end
  end
else
end

nbins = round(width/bin);
if isempty(TS1) || isempty(TS2)
    msgbox('No spikes in one of the chosen clusters.')
else
    [y, x] = CrossCorr(TS1, TS2, bin, nbins);
    bar(x,y,'LineStyle','none')
    xlabel('msec')
    ylabel('Count')
    title('XCorr');
    set(gca,'fontname','arial');
end

if nargout > 0
    Y = y;
    X = x;
else
end

if s == 1
    saveas(gcf, genvarname([n1(1:end-4),'_',n2(1:end-4),'_XCorr']), type);
elseif s ~= 1
end 
    
end