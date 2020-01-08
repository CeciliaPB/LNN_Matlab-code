function OEtoDat(varargin)

% Convert OpenEphys .continuous files to .dat files
%
% INPUTS: 
%   Varargin: 
%   - 'processor': processor number. Default 101.
% 
% EXAMPLES
% OEtoDat;
% OEtoDat('processor',101);
%
% -------------------------------------------------------------------------
% Modified by Cecília Pardo-Bellver, 2019
% Laboratory of Network Neurophysiology
% Institute of Experimantal Medicine, Hungary.
% -------------------------------------------------------------------------
processor = 101;

if nargin
  for i=1:size(varargin,2)
    switch varargin{i}
      case 'processor'
        processor = varargin{i+1};
      otherwise
    end
  end
else
end

common_avg = common_avg_1ch(cd, 'processor',processor, 32);

for ii = 1:32
[data, ts, ~] = load_open_ephys_data([cd '\' num2str(processor) '_CH' num2str(ii) '.continuous']);

data = data - common_avg;

name = ['CH',num2str(ii),'.dat'];
fid = fopen(name,'w+'); 
fwrite(fid,data,'int16'); % Remember the data type to open it.
fclose(fid);

name = 'Time.dat';
fid = fopen(name,'w+'); 
fwrite(fid,ts,'int16'); % Remember the data type to open it.
fclose(fid);

clearvars data ts
end
