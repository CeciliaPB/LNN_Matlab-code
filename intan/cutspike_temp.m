function waveforms = cutspike(data,allspikes,sr,win)
%CUTSPIKE   Extract spike waveforms.

% Waveform window
winp = win * sr;   % waveform window in data points

% Edge effect
allspikes(allspikes<-winp(1)|allspikes>length(data)-winp(2)) = [];  % first and last spikes may be cropped

% Waveform
winx = repmat(allspikes(:)+winp(1),1,sum(abs(winp))) + ...
    repmat(0:sum(abs(winp))-1,length(allspikes),1);
waveforms = nan(size(winx,1),4,size(winx,2));   % waveforms: spikes x channels x time
for iC = 1:4
    cdata = data(:,iC);
    waveforms(:,iC,:) = cdata(winx);   % waveform data
end

winLen = sum(abs(winp));
numSpikes = length(allspikes);
[numDatapoints, numChannels] = size(data);
spkmask = repmat(0:winLen-1,numSpikes,1);   % spike mask
winx = repmat(allspikes(:)+winp(1),1,winLen) + spkmask;
winx = repmat(winx,1,1,numChannels);
winx = permute(winx,[1 3 2]);
zm = zeros(numSpikes,1,winLen);
chmask = cat(2,zm,zm+numDatapoints,zm+2*numDatapoints,zm+3*numDatapoints);
winx = winx + chmask;

waveforms = data(winx);