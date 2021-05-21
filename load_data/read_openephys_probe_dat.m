function read_openephys_probe_dat(probe, varargin)
%READ_OPENEPHYS   Spike detection for open ephys data.
%   READ_OPENEPHYS loads open ephys continuous files, subtracts the average
%   of all channels (referencing), performs filtering, censoring and spike
%   detection of tetrode channels. Files are savet in TT*.mat files for
%   clustering.
%   
%   Probe: 'A1x32', 'Buzs32' 
%
%   Optional input arguments:
%       DATADIR - path name for reading data
%       RESDIR - path name for writing result files
%       TTSPEC - numerical array for subslecting tetrodes
%       RAWDATAFILETAG - if the raw data files have a name tag appended
%
%   Parameter-value input argumets:
%       REFERENCE - allows referenceing options; default: 'common_avg',
%           common average referencing; set to 'none' to apply no offline
%           referencing
%
%   See also OEDISC.
%
% -------------------------------------------------------------------------
%   Balazs Hangya
%   Laboratory of Systems Neuroscience
%   Institute of Experimental Medicine, Budapest, Hungary
%
% Modified by
%   Cecília Pardo-Bellver, 2019
%   Laboratory of Network Neurophysiology
%   Institute of Experimantal Medicine, Hungary.
% -------------------------------------------------------------------------

% Default arguments
prs = inputParser;
addOptional(prs,'datadir',cd,@(s)isempty(s)|isdir(s))   % data directory
addOptional(prs,'resdir','',@(s)isempty(s)|isdir(s))   % results directory
addOptional(prs,'CHspec',32,@isnumeric)   % Number of channels (default: 32 channels)
addOptional(prs,'TTspec',1:8,@isnumeric)   % Number of tetrodes (default, 8 tetrodes)
addOptional(prs,'rawdatafiletag','',@ischar)   % switch for filtering
addOptional(prs,'processor',101,@isnumeric) % Processor number, default 101
addParameter(prs,'reference','common_avg',@(s)ischar(s)|isempty(s))   % switch for referencing
addParameter(prs,'lastTS',0,@isnumeric) %last TimeStamp included (default: end of the recording)
parse(prs,varargin{:})
g = prs.Results;

Th = 35;  % Threshold for the spike detection

info = load_open_ephys_binary([cd '\structure.oebin'],'continuous',1, 'mmap');
switch g.lastTS
    case 0
        lastEL = length(info.Timestamps); % element of the column which belongs to the lastTS
    otherwise
        lastEL = find(info.Timestamps==g.lastTS*30000);
end
switch probe 
    case 'A1x32'
        channel_map = [1 17 16 32 3 19 14 30 9 25 10 20 8 24 2 29 7 26 15 21 11 23 12 28 6 18 13 22 5 27 4 31];
        pairsNeed = 1; % we need all pairs on a linear probe
    case 'Buzs32'
        channel_map = [31 26 27 21 22 23 18 28 29 17 24 32 20 19 25 30 1 2 16 8 3 10 14 9 7 4 15 5 11 13 12 6];
        pairsNeed = 2; % we need only every second pair, because others would be from different shanks
end
        

 % Tetrode organization of the channels
if nargin < 3 || isempty(g.CHspec)
    g.CHspec = 32;
    NumTetrodes = g.CHspec / 4;
    g.TTspec = 1:NumTetrodes;
end

% Directories
if isequal(g.datadir(end),'\')
    g.datadir = g.datadir(1:end-1);
end

SaveFeatures = true;   % save MClust feature files

% Common average reference
switch g.reference
    case 'common_avg'
        common_avg = common_avg_ref_probe_dat(g.datadir,lastEL,g.CHspec,[],g.rawdatafiletag,'processor', g.processor);
    case {'','none'}
        common_avg = 0;
    otherwise
        error('read_openephys: Unknown reference option.')
end

% Spike detection as tetrodes
for iT = g.TTspec
    [tdata, ts] = deal(cell(1,4));
    for iC = 1:4
        Ch1 = (iT - 1) * 4;
        % [tdata{iC}, ts{iC}, info] = load_open_ephys_data([ g.datadir '\' num2str(g.processor) '_CH' num2str(Ch1+iC) g.rawdatafiletag '.continuous']);
        tdata{iC} = info.Data.Data(1).mapped(channel_map(Ch1+iC),1:lastEL);
        tdata{iC} = double(tdata{iC})'*0.195; %intan
        disp (['CH' num2str(channel_map(Ch1+iC)) 'for GR'])
        ts{iC} = double(info.Timestamps)/30000; %in seconds
    end
    data = [tdata{1} tdata{2} tdata{3} tdata{4}];
    data = data - common_avg;
    if isequal(ts{1},ts{2},ts{3},ts{4})   % we assume that continuous channels share timestamps
        tss = ts{1};
    else
        error('Timestamp error.')
    end
    
    [AllTimeStamps, AllWaveForms] = oedisc_probe(data,tss,info.Header.sample_rate,Th,[]);   % filter and detect spikes
    TimeStamps = AllTimeStamps{1};
    WaveForms = AllWaveForms{1};
    GRname = ['GR' num2str(iT)];
    GR_path = fullfile(g.resdir,GRname);
    save(GR_path, 'TimeStamps','WaveForms');
    
    if SaveFeatures   % pre-calculate MClust features
        GRdata.TimeStamps = TimeStamps;
        GRdata.WaveForms = WaveForms;
        openephys_SaveMClustFeatures(GRdata,{'Amplitude';'Energy';'WavePC1';'WavePC2';'Time'},[1 1 1 1],GR_path);
    end
    
    clearvars -except g common_avg filepath resdir SaveFeatures iT Th info channel_map lastEL pairsNeed
end

% Pairs
for iT = 1:pairsNeed:max(g.TTspec) 
        [tdata, ts] = deal(cell(1,4));
    for iC = 1:4
        Ch1 = (iT - 0.5) * 4;
        %[tdata{iC}, ts{iC}, info] = load_open_ephys_data([ g.datadir '\' num2str(g.processor) '_CH' num2str(Ch1+iC) g.rawdatafiletag '.continuous']);
        tdata{iC} = info.Data.Data(1).mapped(channel_map(Ch1+iC),1:lastEL);
        tdata{iC} = double(tdata{iC})'*0.195; %intan
        disp (['CH' num2str(channel_map(Ch1+iC)) 'for PR'])
        ts{iC} = double(info.Timestamps)/30000; %in seconds
    end
    data = [tdata{1} tdata{2} tdata{3} tdata{4}];
    data = data - common_avg;
    if isequal(ts{1},ts{2},ts{3},ts{4})   % we assume that continuous channels share timestamps
        tss = ts{1};
    else
        error('Timestamp error.')
    end
    
    [AllTimeStamps, AllWaveForms] = oedisc_probe(data,tss,info.Header.sample_rate,Th,[]);   % filter and detect spikes
    TimeStamps = AllTimeStamps{1};
    WaveForms = AllWaveForms{1};
    GRname = ['GRinter' num2str(iT)];
    GR_path = fullfile(g.resdir,GRname);
    save(GR_path, 'TimeStamps','WaveForms');
    
    if SaveFeatures   % pre-calculate MClust features
        GRdata.TimeStamps = TimeStamps;
        GRdata.WaveForms = WaveForms;
        openephys_SaveMClustFeatures(GRdata,{'Amplitude';'Energy';'WavePC1';'WavePC2';'Time'},[1 1 1 1],GR_path);
    end
    
    clearvars -except g common_avg filepath resdir SaveFeatures iT Th info channel_map lastEL pairsNeed
    
    
end