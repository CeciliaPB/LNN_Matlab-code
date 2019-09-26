function intandisc(data,thr,resdir)
%INTANDISC   Unit discrimination.
%   INTANDISC(DATA,THR) performs threshold discrimination of continuous
%   unit data (DATA) using the specified threshold (THR). Peak times
%   ('TimeStamps') and spike waveforms ('WaveForms') are saved for each
%   tetrode. A 750 us censored period is applied after each spike. Time
%   window for waveform data is set to -240 to 960 us.
%
%   INTANDISC(DATA,THR,DR) saves the results in the specified directory 
%   (DR). 
%
%   See also READ_INTAN_DATA.

%   Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   9-May-2013

% Results directory
global DATADIR
error(nargchk(2,3,nargin))   % input argument check
if nargin < 3
    resdir = [DATADIR 'Intan' filesep 'nba47' filesep];
end
savestr0 = ['save(''' resdir filesep 'TT'];

% Sampling rate
sr = 25000;   % 25 kHz
deadtime = 0.00075;   % 750 us dead time
dtp = deadtime * sr;   % dead time in data points

% Waveform window
win = [-0.00024 0.00096];   % -240 to 960 us
winp = win * sr;   % waveform window in data points

% Threshold discrimintaion
NumChannels = size(data,2);   % number of channels - 16
NumTetrodes = NumChannels / 4;   % number of tetrodes - 4
tvdisc = cell(1,NumTetrodes);
for iT = 1:NumTetrodes
    [cvdisc tdata] = deal(cell(1,4));
    for iC = 1:4
        chnum = (iT - 1) * 4 + iC;  % current channel
        disp(['tetrode: ' num2str(iT) '   Channel: ' num2str(chnum)])
        cdata = data(:,chnum)';   % data from one channel
        cvdisc{iC} = disc(-cdata,thr);   % discriminate (spike times)
        tdata{iC} = cdata;   % data from the current tetrode
    end
    tvdisc{iT} = cell2mat(cvdisc);
    tvdisc{iT} = sort(tvdisc{iT},'ascend');  % all spike times from one tetrode
        
    % Dead time
    dtv = diff(tvdisc{iT});   % ISI
    censor_inx = dtv < dtp;   % ISI < dead time
    tvdisc{iT}(find(censor_inx)+1) = [];   % remove spikes within the censor period
    
    % Waveform
    winx = repmat(tvdisc{iT}(:)+winp(1),1,sum(abs(winp))) + repmat(0:sum(abs(winp))-1,length(tvdisc{iT}),1);
    wv = nan(size(winx,1),4,size(winx,2));   % waveforms: spikes x channels x time
    for iC = 1:4
        wv(:,iC,:) = tdata{iC}(winx);   % waveform data
    end
    WaveForms = wv; %#ok<NASGU>
    
    % Spike times
    spike_times = tvdisc{iT} / sr;
    TimeStamps = spike_times; %#ok<NASGU>
    savestr = [savestr0 num2str(iT) ''',''WaveForms'',''TimeStamps'');'];
    
    % Save
    eval(savestr)
end