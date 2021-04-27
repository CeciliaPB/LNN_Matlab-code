function read_openephys_CP(varargin)
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
% -------------------------------------------------------------------------
%   Balazs Hangya and Panna Hegedus
%   Laboratory of Systems Neuroscience
%   Institute of Experimental Medicine, Budapest, Hungary
%
%   Modified by Cecília Pardo-Bellver, 2021
%   Laboratory of Network Neurophysiology
%   Institute of Experimantal Medicine, Hungary.
% -------------------------------------------------------------------------

% Default arguments
prs = inputParser;
addOptional(prs,'datadir','n:\Pancsi\HDB12\2017-02-28_17-57-37_HDB12\',@(s)isempty(s)|isdir(s))   % data directory
addOptional(prs,'resdir','',@(s)isempty(s)|isdir(s))   % results directory
addOptional(prs,'TTspec',1:8,@isnumeric)   % selected tetrodes (default: 32 channels, 8 tetrodes)
addOptional(prs,'rawdatafiletag','',@ischar)   % switch for filtering
addParameter(prs,'reference','common_avg',@(s)ischar(s)|isempty(s))   % switch for referencing
parse(prs,varargin{:})
g = prs.Results;

% Tetrodes to convert
NumChannels = 32;   % number of channels
NumTetrodes = NumChannels / 4;   % number of tetrodes
if nargin < 3 || isempty(g.TTspec)
    g.TTspec = 1:NumTetrodes;
end

% Directories
if isequal(g.datadir(end),'\')
    g.datadir = g.datadir(1:end-1);
end

% if isempty(g.resdir)
%     if any(ismember(g.datadir,'-'))
%         sessiontag = 'a';
%         cmps = strsplit(g.datadir,'\');  % path components
%         animalID = cmps{cellfun(@(s)~isempty(s),regexp(cmps,'HDB(\d+)'))};
%         pd = regexp(cmps{end},'(\d+)-(\d+)-(\d+)','tokens');  % extract date
%         sessionID = [pd{1}{1}(3:4) pd{1}{2} pd{1}{3} sessiontag];
%         g.resdir = fullfile(getpref('cellbase','datapath'),animalID,sessionID);
%     else
%         cmps = strsplit(g.datadir,'\');  % path components
%         animalID = cmps{cellfun(@(s)~isempty(s),regexp(cmps,'HDB(\d+)'))};
%         sessionID = cmps{end};
%     end
%     g.resdir = fullfile(getpref('cellbase','datapath'),animalID,sessionID);
%     g.resdir = fullfile('e:\HDBpavlovian_cellbase\',animalID,sessionID);
% end

SaveFeatures = true;   % save MClust feature files

% Common average reference
switch g.reference
    case 'common_avg'
        common_avg = common_avg_ref(g.datadir,32,[],g.rawdatafiletag);
    case {'','none'}
        common_avg = 0;
    otherwise
        error('read_openephys: Unknown reference option.')
end

% Spike detection
for iT = g.TTspec
    [tdata, ts] = deal(cell(1,4));
    for iC = 1:4
        basechannel = (iT - 1) * 4;
        [tdata{iC}, ts{iC}, info] = load_open_ephys_data([ g.datadir '\' '100_CH' num2str(basechannel+iC) g.rawdatafiletag '.continuous']);
    end
    data = [tdata{1} tdata{2} tdata{3} tdata{4}];
    data = data - common_avg;
    if isequal(ts{1},ts{2},ts{3},ts{4})   % we assume that continuous channels share timestamps
        tss = ts{1};
    else
        error('Timestamp error.')
    end
    
    [AllTimeStamps, AllWaveForms] = oedisc(data,tss,info.header.sampleRate,30,[]);   % filter and detect spikes
    TimeStamps = AllTimeStamps{1};
    WaveForms = AllWaveForms{1};
    TTname = ['TT' num2str(iT)];
    TT_path = fullfile(g.resdir,TTname);
    save(TT_path, 'TimeStamps','WaveForms');
    
    if SaveFeatures   % pre-calculate MClust features
        TTdata.TimeStamps = TimeStamps;
        TTdata.WaveForms = WaveForms;
        openephys_SaveMClustFeatures(TTdata,{'Amplitude';'Energy';'WavePC1';'Time'},[1 1 1 1],TT_path)
    end
    
    clearvars -except g common_avg filepath resdir SaveFeatures iT
end