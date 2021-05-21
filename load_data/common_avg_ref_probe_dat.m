function common_avg = common_avg_ref_probe_dat(filepath,lastEL,varargin)
%COMMON_AVG_REF   Common average referencing.
%   AVG = COMMON_AVG_REF(PATH) calculates average of 32 open ephys
%   recording channels for performing offline common average referencing.
%
%   AVG = COMMON_AVG_REF(PATH,CHANNEL_NUMBER,EXCLUDE_CHANNELS) takes
%   optional input arguments for number of recording channels (default, 32)
%   and any potential channels to exclude from averaging.
%
%   See also LOAD_OPEN_EPHYS_DATA.
%
% -------------------------------------------------------------------------
% Modified by Cecília Pardo-Bellver, 2019
% Laboratory of Network Neurophysiology
% Institute of Experimantal Medicine, Hungary.
%
% Based on COMMON_AVG_REF.m
% By Panna Hegedüs, Balazs Hangya 2017/09
% -------------------------------------------------------------------------

% Default arguments
prs = inputParser;
addRequired(prs,'filepath',@(s)exist(s,'dir'))   % data path
addOptional(prs,'channel_number',32,@isnumeric)   % number of recording channels
addOptional(prs,'exclude_channels',[],@isnumeric)   % exclude these channels from averaging
addOptional(prs,'rawdatafiletag','',@ischar)   % tags for raw data filenames
addOptional(prs,'processor',101,@isnumeric) % Processor number, default 101
parse(prs,filepath,varargin{:})
g = prs.Results;

% Calculate sum of all channels
info = load_open_ephys_binary([cd '\structure.oebin'],'continuous',1, 'mmap');
channelsum = [];
for iC = 1:g.channel_number
    if ~ismember(iC,g.exclude_channels)   % exclude channels
        data = info.Data.Data(1).mapped(iC,1:lastEL);
        data = double(data)'*0.195; %intan
        pcs = horzcat(channelsum,data);
        channelsum = sum(pcs,2);
        disp(['CAR - ' num2str(iC)])
    end
end

% Calculate average
NumChannels = g.channel_number - length(g.exclude_channels);
common_avg = channelsum / NumChannels;   % average
common_avg = [common_avg common_avg common_avg common_avg];  % repeat four times to subtract from tetrode data