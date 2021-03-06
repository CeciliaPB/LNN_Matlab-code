function [HISTvals,X] = ACorrelogram(neuron, varargin)

% This function cauculates the autocorrelation of a neuron neurons.
%
% INPUT: 
%   - neuron: .mat file with thw timestamps of the neuron. ex. 'GR1_1.mat'
%   Varargin
%   - 'bins': bin size of choice in msec. Default 1.
%   - 'width': msec on each side of the xcorr plot. Default 250.
%   - 's': saves the figure in the specified file format (jpg, png,
%   svg, etc.). See saveas help for more info. Default NOT save.
%
% OUTPUT:
%   - HISTvals: histogram bars of the plot
%   - X: x dim of the plot, in msec.
%   - without output arguments, the function plots the xcorr
% 
% Examples: 
% [HISTvals, X] = ACorrelogram('GR1_1.mat', 'bin', 1, 'width', 250, 's', 'jpg');
% ACorrelogram('GR1_1.mat', 's', 'jpg');
% 
% -------------------------------------------------------------------------
% Cec�lia Pardo-Bellver, 2019
% Laboratory of Network Neurophysiology
% Institute of Experimental Medicine, Hungary.
%
% Uses code from MClust, MEX-file AutoCorr.
% -------------------------------------------------------------------------

load(neuron,'TS');
TS1 = TS;

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


if isempty(TS1)
    msgbox('No spikes in the chosen clusters.')
else
    [histvals,x] = AutoCorr(TS1, bin, width);
    bar(x,histvals,'LineStyle','none')
    xlabel('ms')
    xlabel(['ms(' num2str(bin) 'ms binsize)']);
    ylabel('rate');
    title(['Autocorrelation - Neuron ',neuron(1:end-4)],'Interpreter','none');
    set(gca,'fontname','arial');
    
if nargout > 0
    HISTvals = histvals;
    X = x;
else
end

if s == 1
    saveas(gcf, genvarname([neuron(1:end-4),'_ACorr']), type);
elseif s ~= 1
end 

end