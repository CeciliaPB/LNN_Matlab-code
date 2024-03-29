function lightpsth_2laser(sessionpath)
%LIGHTPSTH   Spike times aligned to photostimulation.
%   LIGHTPSTH(SESSIONPATH) calculates peri-event time histogram aligned to
%   the onset of photostimulation pulse trains for unclustered tetrode
%   data.
%
%   Saa also TAGGEDPROP.
% -------------------------------------------------------------------------
%   Cec�lia Pardo-Bellver, 2021
%   Laboratory of Network Neurophysiology
%   Institute of Experimental Medicine, Hungary.
%   Based on lightpsth by Balazs Hangya.
% -------------------------------------------------------------------------

% List of TT files
dr = dir;
files = {dr.name};
TTpattern = getpref('cellbase','cell_pattern');   % find tetrode files with a cellbase-defined naming convention
TTinx = regexp(files,[TTpattern '\d\.mat']);
TTinx = cellfun(@(s)~isempty(s),TTinx);
TTfiles = files(TTinx);

% Load photostimulation time stamps
fnm = fullfile(sessionpath,'TTLOnOff.mat');
if ~exist(fnm,'file')
    dr = dir('**/all_channels.events');
    [~, pulseon, ~] = load_events(dr.folder);   % convert PulsePal events
else
    dr = dir('**/all_channels.events');
end
if ~exist('pulseon','var')
    load(fnm,'pulseon');
end

% Load spike times
NumTetrodes = length(TTfiles);
wn = [-20 100];   % PSTH window: -20 to 100 ms 
mwn = max(abs(wn));   % maximal lag
psth = nan(NumTetrodes,2*mwn+1);
legendstring = cell(1,NumTetrodes);

if      ~isstruct(pulseon)
    for iT = 1:NumTetrodes
        fnm = fullfile(sessionpath,TTfiles{iT});
        load(fnm,'TimeStamps')
        
        % Pseudo-trains
        mx = ceil(TimeStamps(end)*1000);
        pse = zeros(1,mx+1000);   % pseudo-event train, ms resolution
        pse(ceil(pulseon*1000)) = 1;
        psu = zeros(1,mx+1000);   % pseudo-spike train
        psu(ceil(TimeStamps*1000)) = 1;
        
        % PSTH
        [lpsth, lags] = xcorr(psu,pse,mwn);
        psth(iT,:) = lpsth;
        legendstring{iT} = ['TT' num2str(iT)];
    end
    
    % Plot
    H = figure;
    plot(lags,psth');
    legend(TTfiles);
    fnm = [sessionpath '\' 'light_psth.jpg'];   % save
    saveas(H,fnm)
    fnm = [sessionpath '\' 'light_psth.fig'];
    saveas(H,fnm)
    pause(0.5)
    close(H)
    
elseif isstruct(pulseon)
    TTLname = fieldnames(pulseon);
    numTTL = size(TTLname,1);
    for nTTL = 1:numTTL
        newpulse = pulseon.(TTLname{nTTL,1});
        for iT = 1:NumTetrodes
            fnm = fullfile(sessionpath,TTfiles{iT});
            load(fnm,'TimeStamps')

            % Pseudo-trains
            mx = ceil(TimeStamps(end)*1000);
            pse = zeros(1,mx+1000);   % pseudo-event train, ms resolution
            pse(ceil(newpulse*1000)) = 1;
            psu = zeros(1,mx+1000);   % pseudo-spike train
            psu(ceil(TimeStamps*1000)) = 1;

            % PSTH
            [lpsth, lags] = xcorr(psu,pse,mwn);
            psth(iT,:) = lpsth;
            legendstring{iT} = ['TT' num2str(iT)];
        end

        % Plot
        H = figure;
        plot(lags,psth');
        legend(TTfiles);
        fnm = [TTLname{nTTL,1} '_light_psth.jpg'];   % save
        saveas(H,fnm)
        pause(0.5)
        close(H)
    end
end

