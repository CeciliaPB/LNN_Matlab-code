function read_nlx(filepath,resdir,TTspec)
%READ_NLX   Spike detection for continuous Neuralynx data.
%   READ_NLX loads Neuralynx continuous files, subtracts the average
%   of all channels (referencing), performs filtering, censoring and spike
%   detection of tetrode channels. Files are savet in TT*.mat files for
%   clustering.
%
%   See also READ_OPENEPHYS and OEDISC.

%   Balazs Hangya and Panna Hegedus
%   Laboratory of Systems Neuroscience
%   Institute of Experimental Medicine, Budapest, Hungary

% Tetrodes to convert
NumChannels = 32;   % number of channels
NumTetrodes = NumChannels / 4;   % number of tetrodes
if nargin < 3 || isempty(TTspec)
    TTspec = 1:NumTetrodes;
end

% Directories
if nargin < 1
    filepath = 'G:\_data\_recordings\HDB23\2017-12-05_14-08-35_HDB23';
    filepath = 'n:\Pancsi\HDB12\2017-02-28_17-57-37_HDB12\';
end
if isequal(filepath(end),'\')
    filepath = filepath(1:end-1);
end

if ~exist('resdir','var') || isempty(resdir)
    if any(ismember(filepath,'-'))
        sessiontag = 'a';
        cmps = strsplit(filepath,'\');  % path components
        animalID = cmps{cellfun(@(s)~isempty(s),regexp(cmps,'HDB(\d+)'))};
        pd = regexp(cmps{end},'(\d+)-(\d+)-(\d+)','tokens');  % extract date
        sessionID = [pd{1}{1}(3:4) pd{1}{2} pd{1}{3} sessiontag];
        resdir = fullfile(getpref('cellbase','datapath'),animalID,sessionID);
    else
        cmps = strsplit(filepath,'\');  % path components
        animalID = cmps{cellfun(@(s)~isempty(s),regexp(cmps,'HDB(\d+)'))};
        sessionID = cmps{end};
    end
    resdir = fullfile(getpref('cellbase','datapath'),animalID,sessionID);
    resdir = fullfile('e:\HDBpavlovian_cellbase\',animalID,sessionID);
end

SaveFeatures = true;   % save MClust feature files

% Common average reference
common_avg = common_avg_ref(filepath,32);

% Spike detection
for iT = TTspec
    [tdata, ts] = deal(cell(1,4));
    for iC = 1:4
        basechannel = (iT - 1) * 4;
        [tdata{iC}, ts{iC}, sr] = load_nlx_data(fullfile(filepath,['CSC' num2str(basechannel+iC) '.ncs']));
    end
    data = [tdata{1} tdata{2} tdata{3} tdata{4}];
    data = data - common_avg - eps;
    if isequal(ts{1},ts{2},ts{3},ts{4})   % we assume that continuous channels share timestamps
        tss = ts{1};
    else
        error('Timestamp error.')
    end
    [AllTimeStamps, AllWaveForms] = oedisc(data,tss,sr,30,[]);   % filter and detect spikes
    TimeStamps = AllTimeStamps{1};
    WaveForms = AllWaveForms{1};
    TTname = ['TT' num2str(iT)];
    TT_path = fullfile(resdir,TTname);
    save(TT_path, 'TimeStamps','WaveForms');
    
    if SaveFeatures   % pre-calculate MClust features
        TTdata.TimeStamps = TimeStamps;
        TTdata.WaveForms = WaveForms;
       % openephys_SaveMClustFeatures(TTdata,{'Amplitude';'Energy';'WavePC1';'Time'},[1 1 1 1],TT_path)
       openephys_SaveMClustFeatures(TTdata,{'Amplitude';'Energy';'Time'},[1 1 1 1],TT_path)
    end
    
    clearvars -except common_avg filepath resdir SaveFeatures iT rawdatafiletag
end

% -------------------------------------------------------------------------
function [data, time, sr] = load_nlx_data(cscname)

% Load continuous data
disp(['Reading ' cscname])
[TimeStamps_CSC, ~, ~, ~, Samples, NlxHeader] = ...
    Nlx2MatCSC(cscname,[1 1 1 1 1],1,1,1);
TimeStamps_CSC = TimeStamps_CSC / 1e6;   % convert time to seconds
data = -Samples(1:end)';

NHs = NlxHeader{cellfun(@(s)~isempty(s),strfind(NlxHeader,'ADBitVolts'))};   % header scale variable
NHsv = regexp(NHs,'\d\.\d*','match');
ADBitVolts = str2double(NHsv{1});
data = data * ADBitVolts * 1e06;
sr = 512 / (mean(diff(TimeStamps_CSC)));   % sampling rate
dt = 1 /sr;  % sampling time
time = TimeStamps_CSC(1):dt:TimeStamps_CSC(1)+dt*(length(data)-1);  % time vector (even)
starttime = TimeStamps_CSC(1);   % first time stamp

% -------------------------------------------------------------------------
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

% Panna Hegedüs, Balazs Hangya 2017/09

% Default arguments
prs = inputParser;
addRequired(prs,'filepath',@(s)exist(s,'dir'))   % data path
addOptional(prs,'channel_number',32,@isnumeric)   % number of recording channels
addOptional(prs,'exclude_channels',[],@isnumeric)   % exclude these channels from averaging
parse(prs,filepath,varargin{:})
g = prs.Results;

% Calculate sum of all channels
channelsum = [];
for iC = 1:g.channel_number
    if ~ismember(iC,g.exclude_channels)   % exclude channels
        data = load_nlx_data(fullfile(filepath,['CSC' num2str(iC) '.ncs']));
        pcs = horzcat(channelsum,data);
        channelsum = sum(pcs,2);
    end
end

% Calculate average
NumChannels = g.channel_number - length(g.exclude_channels);
common_avg = channelsum / NumChannels;   % average
common_avg = [common_avg common_avg common_avg common_avg];  % repeat four times to subtract from tetrode data