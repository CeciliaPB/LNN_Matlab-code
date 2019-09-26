function [AllTimeStamps, AllWaveForms] = oedisc_1ch(data,ts,sr,thr,varargin)
%OEDISC   Unit discrimination.
%   [T W] = OEDISC(DATA,TS,SR,THR) performs threshold discrimination of
%   continuous timestamped (TS) unit data (DATA) sampled at the rate SR
%   using the specified threshold (THR). Peak times (T, 'TimeStamps') and
%   spike waveforms (W, 'WaveForms') are saved for each tetrode. A 750 us
%   censored period is applied and the larger spike is kept. Time window
%   for waveform data is set to -300 to 1000 us.
%
%   OEDISC(DATA,TS,SR,THR,DR) saves the results in the specified directory 
%   (DR). If DR is empty, the data is not saved.
%
%   See also READ_INTAN_DATA.

%   Balazs Hangya, Institute of Experimental Medicine
%   hangya.balazs@koki.mta.hu
%   30-Jan-2018

% Default arguments
prs = inputParser;
addRequired(prs,'data',@isnumeric)   % raw or filtered data
addRequired(prs,'ts',@isnumeric)   % timestamps
addRequired(prs,'sr',@isnumeric)   % sampling rate
addRequired(prs,'thr',@isnumeric)   % discrimination threshold
addOptional(prs,'resdir',['C:' filesep 'Balazs' filesep '_data' filesep 'Intan' filesep 'Tina' filesep],...
    @(s)ischar(s)|isempty(s))   % results directory
addParameter(prs,'Filtering','enable',@(s)ischar(s)|...
    ismember(s,{'disable','enable'}))   % switch for filtering
parse(prs,data,ts,sr,thr,varargin{:})
g = prs.Results;

% File name
savestr0 = ['save(''' g.resdir filesep 'TT'];

% Sampling rate
nqf = sr / 2;   % Nyquist freq.
deadtime = 0.00075;   % 750 us dead time
dtp = deadtime * sr;   % dead time in data points

% Waveform window
win = [-0.0003 0.0006];   % -300 to 600 us
winp = round(win*sr);   % waveform window in data points

% Threshold discrimintaion
NumChannels = size(data,2);   % number of channels - 16
[tvdisc, tpeaks, AllTimeStamps, AllWaveForms] = deal(cell(1,NumChannels));
for iT = 1:NumChannels
    [cvdisc, tdata, cpeaks] = deal(cell(1,4));
    
    chnum = iT ;  % current channel
    disp(['tetrode: ' num2str(iT) '   Channel: ' num2str(chnum)])
    cdata = data(:,chnum)';   % data from one channel
    if isequal(g.Filtering,'enable')
        [b,a] = butter(3,[700 7000]/nqf,'bandpass');   % Butterworth filter
        unit = filter(b,a,cdata);  % filter
    elseif isequal(g.Filtering,'disable')
        unit = cdata;
    else
        error('oedisc:InputArg','Unsupported input argument for filtering.')
    end
    [cvdisc{iT}, cpeaks{iT}] = disc(-unit,thr);   % discriminate (spike times)
    tdata{iT} = -unit;   % data from the current tetrode
   
    tvdisc{iT} = cell2mat(cvdisc);
    tpeaks{iT} = cell2mat(cpeaks);
    [tvdisc{iT}, inx] = sort(tvdisc{iT},'ascend');  % all spike times from one tetrode
    tpeaks{iT} = tpeaks{iT}(inx);
        
    % Dead time
    dtv = diff(tvdisc{iT});   % ISI
    dtpk = diff(tpeaks{iT});   % comparison of neighboring peaks
    while any(dtv<dtp)
        censor_inx = dtv < dtp;   % ISI < dead time
        peak_comp = dtpk > 0;
        delete_inx1 = censor_inx & peak_comp;
        delete_inx2 = [0 censor_inx & ~peak_comp];
        tvdisc{iT}([find(delete_inx1) find(delete_inx2)]) = [];
        tpeaks{iT}([find(delete_inx1) find(delete_inx2)]) = [];
        dtv = diff(tvdisc{iT});   % ISI
        dtpk = diff(tpeaks{iT});   % comparison of neighboring peaks
    end
    tvdisc{iT}(tvdisc{iT}<=-winp(1)|tvdisc{iT}>=size(data,1)-winp(2)) = [];   % we may lose some spikes near the ends of the file
    
    % Waveform
    winx = repmat(tvdisc{iT}(:)+winp(1),1,sum(abs(winp))) + repmat(0:sum(abs(winp))-1,length(tvdisc{iT}),1);
    wv = nan(size(winx,1),4,size(winx,2));   % waveforms: spikes x channels x time
    wv(:,iT,:) = tdata{iT}(winx);   % waveform data
    WaveForms = wv;
    AllWaveForms{iT} = WaveForms;
    
    % Spike times
    spike_times = ts(tvdisc{iT});
    TimeStamps = spike_times;
    AllTimeStamps{iT} = TimeStamps;
    savestr = [savestr0 num2str(iT) ''',''WaveForms'',''TimeStamps'');'];
    
    % Save
    if ~isempty(g.resdir)
        eval(savestr)
    end
end