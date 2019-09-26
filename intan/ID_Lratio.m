function [DM LrC] = ID_Lratio(spiketimes,features,valid_channels)
%LRATIO2   L-ratio and Isolation Distance.
%   [DM LRC] = LRATIO2(CELLID,'FEATURE_NAMES',FN,'VALID_CHANNELS',VC)
%   calculates cluster quality measures (Isolation Distance, DM;
%   L-ratio,LRC) for a given cell (CELLID) using the given clustering
%   features (FN; default: {'Energy' 'WavePC1'}). VC is a 4-element 0-1
%   array determining which of the four tetrode channels were functional.
%   If it is omitted, the program defines it based on valid channels having
%   nonzero waveform energy.
%
%   [DM LRC VALID_CHANNELS] = LRATIO2(CELLID,'FEATURE_NAMES',FN,'VALID_CHANNELS',VC)
%   also returns channel validity.
%
%   See also MAHAL.

%   Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com

% Channel validity
if nargin < 3
    valid_channels = [1 2 3 4];
end

% Isolation Distance
XC = X(:,inx);
D = mahal(X',XC');
DC = D(cinx);
sDC = sort(DC);
linx = length(inx);
if linx <= length(sDC)
    DM = sDC(length(inx));
else
    DM = Inf;   % more spikes in the cluster than outside of it
end

% L-ratio
if ~isempty(DC)
    df = size(X,1);
    LC = sum(1-chi2cdf(DC,df));
    nc = size(XC,2);
    LrC = LC / nc;
else
    LrC = Inf;  % no spikes outside of the cluster: multiunit
end

% -------------------------------------------------------------------------
function calculate_features(sessionpath,propfn_path,feature_names,basename,valid_channels)

% Create MClust variables
global MClust_FDdn
global MClust_ChannelValidity
global MClust_NeuralLoadingFunction
global MClust_TText
global MClust_FDext
global MClust_TTdn
global MClust_TTfn

MClust_FDext = '.fd';
MClust_TText = '.ntt';
MClust_TTfn = basename;
[t1, t2] = strtok(fliplr(which('MClust')),filesep);
MClust_Directory = fliplr(t2);
MClust_FDdn = propfn_path;
MClust_TTdn = sessionpath;
MClust_ChannelValidity = valid_channels;
MClust_NeuralLoadingFunction = char([MClust_Directory 'LoadingEngines\LoadTT_NeuralynxNT']);

% Calculate features
CalculateFeatures(basename,feature_names)