function waveforms = cutspike(data,allspikes,sr,win)
%CUTSPIKE   Extract spike waveforms.
%   WAVEFORMS = CUTSPIKE(DATA,ALLSPIKES,SR,WIN) extracts the spike
%   waveforms around the time stamps in ALLSPIKES from the DATA (sampled at
%   SR sampling rate). DATA should be in samples x channels format. Window
%   size is given by WIN in seconds. WAVEFORMS are returned in spikes x
%   channels x time format.
%
%   Note that as CUTSPIKE was meant to be a low level utility function,
%   ALLSPIKES should be given in data points (not converted into seconds).
%
%   See also PEAKDETECT and CENSOR.

%   Balazs Hangya, TENSS 2016
%   hangya.balazs@koki.mta.hu

% Waveform window
winp = win * sr;   % waveform window in data points

% Edge effect
allspikes(allspikes<-winp(1)|allspikes>length(data)-winp(2)) = [];  % first and last spikes may be cropped

% Size variables
winLen = sum(abs(winp));   % number of time points in a spike
numSpikes = length(allspikes);   % number of spikes
[numDatapoints, numChannels] = size(data);   % total number of data points per channel, number of channels
if numChannels > 128
    warning('If you have less then 128 channels, your data matrix is likely transposed.')
end

% Waveforms
spkmask = repmat(0:winLen-1,numSpikes,1);   % spike mask
winx = repmat(allspikes(:)+winp(1),1,winLen) + spkmask;   % spike windows
winx = repmat(winx,1,1,numChannels);
winx = permute(winx,[1 3 2]);
chmask = repmat((0:numChannels-1)*numDatapoints,numSpikes,1,winLen);
winx = winx + chmask;   % spike windows for all channels

waveforms = data(winx);

% This is how it would look with a channel loop
% winx = repmat(allspikes(:)+winp(1),1,sum(abs(winp))) + ...
%     repmat(0:sum(abs(winp))-1,length(allspikes),1);
% waveforms = nan(size(winx,1),4,size(winx,2));   % waveforms: spikes x channels x time
% for iC = 1:numChannels
%     cdata = data(:,iC);
%     waveforms(:,iC,:) = cdata(winx);   % waveform data
% end