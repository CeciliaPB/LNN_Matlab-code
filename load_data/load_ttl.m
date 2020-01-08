function TTL = load_ttl(varargin)

% This function reads a file (.continuous from OpenEphys) or data+time
% variables (.mat generated from an .rhd from Intan) or data+time in the
% Workspace (e. g. audio file) and detects values above a threshold as
% timestamps for the TTL. In the 
% 
% INPUTS: 
%   Varargin: 
%   - file: signal to analyse, string. e. g. '101_ADC1.continuous'
%   - 'origin': where the  signal comes from. CASE 1: origin = 'OE',
%   recorded with OpenEphys; CASE 2: origin = 'IN', recorded with Intan
%   Board, channels saved from the .rhd file; CASE 3: origin = 'WS', for
%   data in the Workspace. Default 'OE'.
%   - 'n': number of TTLs 
%   - 'mPP': MinPeakProminence, theshold for the detection. Default 0.5.
%   The data is normalised in a 0 to 1 range. If you are not sure about the
%   threshold: 1) read the data to Workspace; 2) normalize it: data2 =
%   normalize(data,'range',[0 1]); 3) plot(data2); 4) Choose value.
%   - 'mPD': MinPeakDistance, min time between peaks, in seconds. Default
%   0.1 
%   - 'mPW': MinPeakWidth, min width of the peak, in seconds. Default 0.01 
%   - 's': The TTL is automatically saved in the current folder. Default
%   not save.
%
% OUTPUTS: 
%   - TTL: a vector list of the timestamps of the TTLs.
% 
% Examples: 
% TTL = load_ttl('101_ADC1.continuous');
% TTL = load_ttl('TTL.mat','origin','IN','n',60,'mPP',0.5,...
%       'mPD',0.1,'s');
% TTL = load_ttl('origin','WS','n',60,'mPP',0.5,...
%       'mPD',0.1,'s');
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
mPW = 0.01;
s = 0;
origin = 'OE';

% Params, read from varargin
if nargin
  for i=1:size(varargin,2)
    switch varargin{i}
      case 'origin'
        origin = varargin{i+1};
      case 'n'
        n = varargin{i+1};
      case 'mPP'
        mPP = varargin{i+1};
      case 'mPD'
        mPD = varargin{i+1};
      case 'mPW'
        mPW = varargin{i+1};
      case 's'
        s = 1;
      otherwise
        file = varargin{i};
    end
  end
else
end

switch origin
    case 'OE'
        [data, timestamps, ~] = load_open_ephys_data(file); % Load .continuous files
        data2 = normalize(data,'range',[0 1]);
        data2 = data2>mPP;

        [~,locs] = findpeaks(double(data2),timestamps,'MinPeakProminence'...
            ,mPP,'MinPeakDistance',mPD,'MinPeakWidth',mPW,'NPeaks',n);
        TTL = locs;

        if s == 1
        save('TTL.mat','TTL');
        elseif s ~= 1
        end
        
    case 'IN'
        TTL_id = input('TTL? ENTER NAME (of the variable containing the TTL):','s');
        data = load(file,TTL_id);
        feval(@()assignin('caller','data',  data.(TTL_id)));
        t_dig = 't_dig';
        timestamps = load(file,'t_dig');
        feval(@()assignin('caller','timestamps',  timestamps.(t_dig)));
        data2 = normalize(data,'range',[0 1]);
        data2 = data2>mPP;
        
        [~,locs] = findpeaks(double(data2),timestamps,'MinPeakProminence',mPP,...
            'MinPeakDistance',mPD,'MinPeakWidth',mPW,'NPeaks',n);
        TTL = locs;
        
        if s == 1
        save('TTL.mat','TTL');
        elseif s ~= 1
        end
        
    case 'WS'
        data = input('VARIABLE TO USE AS DATA:'); % e. g. If you tipe audio, 
        % it looks for a variable in the WS with this name to use as your data.
        timestamps = input('VARIABLE TO USE AS TIME:'); %e. g. If you tipe time, 
        % it looks for a variable in the WS with this name to use as your time.
        data2 = normalize(data,'range',[0 1]);
        data2 = data2>mPP;

        [~,locs] = findpeaks(double(data2),timestamps,'MinPeakProminence'...
            ,mPP,'MinPeakDistance',mPD,'MinPeakWidth',mPW,'NPeaks',n);
        TTL = locs;

        if s == 1
        save('TTL.mat','TTL');
        elseif s ~= 1
        end

end
end
