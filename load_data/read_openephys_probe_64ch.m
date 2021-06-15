function read_openephys_probe_64ch(varargin)
%READ_OPENEPHYS   Spike detection for open ephys data.
%   READ_OPENEPHYS loads open ephys continuous files, subtracts the average
%   of all channels (referencing), performs filtering, censoring and spike
%   detection of tetrode channels. Files are savet in TT*.mat files for
%   clustering.
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
%   Balazs Hangya
%   Laboratory of Systems Neuroscience
%   Institute of Experimental Medicine, Budapest, Hungary
%
% -------------------------------------------------------------------------
% To use with probe recordings and Add Option to change the number of the
% processor.
% Modified by Cecília Pardo-Bellver, 2019
% Laboratory of Network Neurophysiology
% Institute of Experimental Medicine, Hungary.
% -------------------------------------------------------------------------

% Default arguments
prs = inputParser;
addOptional(prs,'datadir',cd,@(s)isempty(s)|isdir(s))   % data directory
addOptional(prs,'resdir','',@(s)isempty(s)|isdir(s))   % results directory
addOptional(prs,'CHspec',64,@isnumeric)   % Number of channels (default: 32 channels)
addOptional(prs,'TTspec',1:16,@isnumeric)   % Number of tetrodes (default, 8 tetrodes)
addOptional(prs,'rawdatafiletag','',@ischar)   % switch for filtering
addOptional(prs,'processor',101,@isnumeric) % Processor number, default 101
addParameter(prs,'reference','common_avg',@(s)ischar(s)|isempty(s))   % switch for referencing
parse(prs,varargin{:})
g = prs.Results;

Th = 35;  % Threshold for the spike detection

% Tetrode organization of the channels
if nargin < 3 || isempty(g.CHspec)
    g.CHspec = 64;
    NumTetrodes = g.CHspec / 4;
    g.TTspec = 1:NumTetrodes;
end
NumTetrodes = g.CHspec / 4;
g.TTspec = 1:NumTetrodes;

% Directories
if isequal(g.datadir(end),'\')
    g.datadir = g.datadir(1:end-1);
end

SaveFeatures = true;   % save MClust feature files

% Common average reference
switch g.reference
    case 'common_avg'
        common_avg = common_avg_ref_probe(g.datadir,g.CHspec,[],g.rawdatafiletag,'processor', g.processor);
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
        [tdata{iC}, ts{iC}, info] = load_open_ephys_data([ g.datadir '\' num2str(g.processor) '_CH' num2str(Ch1+iC) g.rawdatafiletag '.continuous']);
    end
    data = [tdata{1} tdata{2} tdata{3} tdata{4}];
    data = data - common_avg;
    if isequal(ts{1},ts{2},ts{3},ts{4})   % we assume that continuous channels share timestamps
        tss = ts{1};
    else
        error('Timestamp error.')
    end
    
    [AllTimeStamps, AllWaveForms] = oedisc_probe(data,tss,info.header.sampleRate,Th,[]);   % filter and detect spikes
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
    
    clearvars -except g common_avg filepath resdir SaveFeatures iT Th
end

% Pairs
for iT = 1:max(g.TTspec)-1
    
    % Channel 1
    Ch1 = (iT * 4) ;
    [tdata1, ts1, ~] = load_open_ephys_data([ g.datadir '\' num2str(g.processor) '_CH' num2str(Ch1) g.rawdatafiletag '.continuous']);
    
    % Channel 2
    Ch2 = Ch1 + 1 ;
    [tdata2, ts2, info] = load_open_ephys_data([ g.datadir '\' num2str(g.processor) '_CH' num2str(Ch2) g.rawdatafiletag '.continuous']);
    
    data = [tdata1 tdata2];
    data = data - common_avg(:,1:2);
    if isequal(ts1,ts2)   % we assume that continuous channels share timestamps
        tss = ts1;
    else
        error('Timestamp error.')
    end
    
    [AllTimeStamps, AllWaveForms] = oedisc_probe(data,tss,info.header.sampleRate,Th,[]);   % filter and detect spikes
    TimeStamps = AllTimeStamps{1};
    WaveForms = AllWaveForms{1};
    GRname = ['PR' num2str(iT+16)];
    GR_path = fullfile(g.resdir,GRname);
    save(GR_path, 'TimeStamps','WaveForms');
    
    if SaveFeatures   % pre-calculate MClust features
        GRdata.TimeStamps = TimeStamps;
        GRdata.WaveForms = WaveForms;
        openephys_SaveMClustFeatures(GRdata,{'Amplitude';'Energy';'WavePC1';'WavePC2';'Time'},[1 1 0 0],GR_path);
    end
    
    clearvars -except g common_avg filepath resdir SaveFeatures iT Th
end