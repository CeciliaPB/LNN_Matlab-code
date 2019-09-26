function read_openephys_dat(varargin)
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

%   Balazs Hangya
%   Laboratory of Systems Neuroscience
%   Institute of Experimental Medicine, Budapest, Hungary

% Default arguments
prs = inputParser;
addOptional(prs,'datadir','n:\Pancsi\HDB12\2017-02-28_17-57-37_HDB12\',@(s)isempty(s)|isdir(s))   % data directory
addOptional(prs,'resdir','',@(s)isempty(s)|isdir(s))   % results directory
addOptional(prs,'TTspec',1:8,@isnumeric)   % selected tetrodes (default: 32 channels, 8 tetrodes)
addOptional(prs,'rawdatafiletag','',@ischar)   % switch for filtering
addParameter(prs,'reference','common_avg',@(s)ischar(s)|isempty(s))   % switch for referencing
parse(prs,varargin{:})
g = prs.Results;

% Tetrode organization of the channels
NumChannels = 32;   % number of channels
NumTetrodes = NumChannels / 4;   % number of tetrodes
if nargin < 3 || isempty(g.TTspec)
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
    
    path = pwd;
    cd(g.resdir)
    
    dat = ['T',num2str(iT),'_data.dat'];
    %initialize a new .dat file
    fid = fopen(dat,'w+');    
    % write data to raw .dat file
    fwrite(fid,data(:,:),'int16');
    % close dat file
    fclose(fid);
    
    ts = ['T','_ts.dat'];
    %initialize a new .dat file
    fid = fopen(ts,'w+');    
    % write data to raw .dat file
    fwrite(fid,tss(:,:),'int16');
    % close dat file
    fclose(fid);
    
    cd(path)
    
    clearvars -except g common_avg filepath resdir SaveFeatures iT
end