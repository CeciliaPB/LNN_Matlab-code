function common_avg = common_avg_ref(filepath,varargin)
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
%   Panna Heged�s, Balazs Hangya 2017/09
%
%   Modified by Cec�lia Pardo-Bellver, 2021
%   Laboratory of Network Neurophysiology
%   Institute of Experimental Medicine, Hungary.
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
channelsum = [];
for iC = 1:g.channel_number
    if ~ismember(iC,g.exclude_channels)   % exclude channels
        data = load_open_ephys_data([filepath '\' num2str(g.processor) '_CH' num2str(iC) g.rawdatafiletag '.continuous']);
        pcs = horzcat(channelsum,data);
        channelsum = sum(pcs,2);
    end
end

% Calculate average
NumChannels = g.channel_number - length(g.exclude_channels);
common_avg = channelsum / NumChannels;   % average
common_avg = [common_avg common_avg common_avg common_avg];  % repeat four times to subtract from tetrode data