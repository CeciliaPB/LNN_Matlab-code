function TTL = load_ttl(file,varargin)

% This function reads a .continuous file or data+time variables in the
% Workspace and detects values above a threshold as timestamps for the TTL
% 
% INPUTS: 
%   - file: signal to analyse, string. Varargins - n: number of TTLs 
%   Varargin:
%   - 'mPP': MinPeakProminence, theshold for the detection. Default 0.5. The
%   data is normalised in a -1 to 1 range. If you are not sure about the
%   threshold: 1) read the data to Workspace, 2) normalize it 'data2 =
%   normalize(data,'range',[-1 1])', 3) plot(data2), 4) Choose value.
%   - 'mPD': MinPeakDistance, min time between peaks, in seconds. Default 0.1 
%   - 's': The TTL is automatically saved in the current folder
%
% OUTPUTS: 
%   - TTL: a vector list of the timestamps of the TTLs.
% 
% Examples: 
% TTL = load_ttl('101_ADC1.continuous');
% TTL = load_ttl('101_ADC1.continuous','n',60,'mPP',0.5,'mPD',0.1,'s');
% 
% -------------------------------------------------------------------------
% Cecília Pardo-Bellver, 2019
% Laboratory of Network Neurophysiology
% Institute of Experimantal Medicine, Hungary.
%
% MATLAB toolbox: Signal Processing Toolbox.
% -------------------------------------------------------------------------

% Params, defalut values:
n = 1000;
mPP = 0.5;     
mPD = 0.1;
s = 0;

% Params, read from varargin
if nargin
  for i=1:size(varargin,2)
    switch varargin{i}
      case 'n'
        n = varargin{i+1};
      case 'mPP'
        mPP = varargin{i+1};
      case 'mPD'
        mPD = varargin{i+1};
      case 's'
        s = 1;
      otherwise
        
    end
  end
else
end

[data, timestamps, ~] = load_open_ephys_data(file); % Load .continuous files
% % Uncomment to use variable in the Workspace as data. 
% data = audio; % audio, is your data variable name;
% timestamps = t; % t, is your time variable name;
data2 = normalize(data,'range',[-1 1]);
data2 = data2>mPP;

[~,locs] = findpeaks(double(data2),timestamps,'MinPeakProminence',mPP,'MinPeakDistance',mPD,'NPeaks',n);

TTL = locs;

if s == 1
save('TTL.mat','TTL');
elseif s ~= 1
end

end