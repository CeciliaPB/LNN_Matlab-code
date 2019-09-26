function [PeakData, PeakNames,PeakPars] = feature_Amplitude(V, ttChannelValidity, Params)

% MClust
% [PeakData, PeakNames] = feature_Amplitude(V, ttChannelValidity)
% Calculate peak feature max value for each channel
%
% INPUTS
%    V = TT tsd
%    ttChannelValidity = nCh x 1 of booleans
%
% OUTPUTS
%    Data - nSpikes x nCh peak values
%    Names - "Peak6to11: Ch"
%

TTData = Data(V);

[nSpikes, nCh, nSamp] = size(TTData);


f = find(ttChannelValidity);


PeakData = zeros(nSpikes, length(f));

PeakNames = cell(length(f), 1);
PeakPars = {};
PeakData = squeeze(max(TTData(:, f, :), [], 3)) - squeeze(min(TTData(:, f, :), [], 3));

for iCh = 1:length(f)
   PeakNames{iCh} = ['Amplitude: ' num2str(f(iCh))];
end
