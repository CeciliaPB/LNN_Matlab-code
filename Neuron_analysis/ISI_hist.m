function [ISIh, ISIbins] = ISI_hist(neuron,varargin)

% This function calculates the interspike interval of the timestamps of the
% neuron.
% 
% INPUTS: 
%   - neuron: .mat file with the timestamps of the neuron. ex. 'GR1_1.mat'
%   Varargin
%   - 'bins': bins used. Afects the size of the bars in the histogram.
%   Default 500.
%   - 's': saves the figure in the specified file format (jpg, png,
%   svg, etc.). See saveas help for more info. Default NOT save.
%
% OUTPUTS: 
%   - ISIh: array of histogram of the neuron timestamps. 
%   - ISIbins: num of bins, in the same units as the timestamps.
%   - without output arguments, the function plots the ISI.
% 
% Examples: 
% [ISIh, ISIbins] = ISI_hist(neuron,'bins',500,'s','jpg')
% ISI_hist(neuron,'bins',500) 
% 
% -------------------------------------------------------------------------
% Cecília Pardo-Bellver, 2019
% Laboratory of Network Neurophysiology
% Institute of Experimental Medicine, Hungary.
%
% Uses code from MClust, MEX-file ndhist.
% -------------------------------------------------------------------------

if ~exist('TS','var') == 1
    load(neuron,'TS');
end

% Default params
epsilon = 1e-100;
bins = 500;
s = 0;

% Params, read from varargin
if nargin
  for i=1:size(varargin,2)
    switch varargin{i}
      case 'bins'
        bins = varargin{i+1};
      case 's'
        s = 1;
        type = varargin{i+1};
      otherwise  
    end
  end
else
end

ISI = diff((TS)/10) + epsilon;
maxLogISI = max(real(log10(ISI)))+1;
minLogISI = floor(min(real(log10(ISI))));

isiH = ndhist(log10(ISI)', bins, minLogISI, maxLogISI);
isiB = logspace(minLogISI,maxLogISI,bins);

% Plot 
plot(isiB, isiH);
set(gca, 'XScale', 'log', 'XLim', [10^minLogISI 10^maxLogISI]);
set(gca, 'YTick', max(isiH));
hold on
plot([1 1], get(gca, 'YLim'), 'r:','LineWidth',2);
title(['InterSpike Interval - Neuron ',neuron(1:end-4)],'Interpreter','none');
hold off

if nargout > 0
    ISIh = isiH;
    ISIbins = isiB;
else
end

if s == 1
    saveas(gcf, genvarname([neuron(1:end-4),'_ISI']), type);
elseif s ~= 1
end    

end
