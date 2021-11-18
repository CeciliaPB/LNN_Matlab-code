function MakeTrialEvents_gambling(sessionpath,varargin)

%MAKETRIALEVENTS2_GONOGO   Synchronize trial events to recording times. 
%	MAKETRIALEVENTS2_GONOGO(SESSIONPATH) loads Neuralynx events and adjusts
%	trial event times (trial-based behavioral data, see
%	SOLO2TRIALEVENTS4_AUDITORY_GONOGO) to the recorded time stamps. This
%	way the neural recordings and behavioral time stamps are in register.
%	Stimulus time TTL pulses are used for synchronization. The synchronized
%	trial events structure is saved under the name 'TrialEvents.mat'. This
%	file becomes the primary store of behavioral data for a particular
%	session; it is retrieved by LOADCB via CELLID2FNAMES. This default
%	file name is one of the preference settings of CellBase - type
%   getpref('cellbase','session_filename').
%
%   MAKETRIALEVENTS2_GONOGO(SESSIONPATH,'StimNttl',TTL) specifies the TTL
%   channel which serves as the basis for synchronization.
%
%   See also LOADCB.
% -------------------------------------------------------------------------
% Cecília Pardo-Bellver, 2021
% Laboratory of Network Neurophysiology
% Institute of Experimental Medicine, Hungary.
% -------------------------------------------------------------------------

eventtype = 'event';
if nargin == 2
    eventtype = varargin{1};
end

% Check session path ------------------------------------------------------
if ~isfolder(sessionpath)
    error('Session path is not valid.');
end

switch  eventtype
    case {'behav','event'}
        % Load Open Ephys events file -------------------------------------
        try
            load([sessionpath '\' 'TTLs.mat'],'TTL');  
        catch
            A = dir('BA*');
            A = A.name;
            [TTL, ~, ~] = load_events([sessionpath '\' A]);
        end
        
        % % Manually choose the TTL ---------------------------------------
        % TTLfileds = fieldnames(event);
        % list = TTLfileds;
        % Newevents = listdlg('PromptString','Choose Behaviour Channel',...
        %     'SelectionMode','single', 'ListString',list);
        % switch Newevents
        %     case 1
        %         event = event.(list{1});
        %     case 2
        %         event = event.(list{2});
        %     case 3
        %         event = event.(list{3});
        %     case 4
        %         event = event.(list{4});
        %     case 5
        %         event = event.(list{5});
        %     case 6
        %         event = event.(list{6});
        %     case 7
        %         event = event.(list{7});
        %     case 8
        %         event = event.(list{8});
        % end
        % TTL = event;
        TTL = TTL.CH8; % Events in CH8 are the behaviour
        
        % Load Trial Events structure -------------------------------------
        SE_filename  = [sessionpath filesep 'TE_behaviour.mat'];
        TE_behaviour = load(SE_filename);
        
        % Synchronisation -------------------------------------------------
        TE2  = TE_behaviour;
        son2 = TTL';   % stimulus onset time recorded by the recording system (OE)
        ts   = TE_behaviour.StimulusOn + TE_behaviour.TrialStart;   % stimulus onset in absolut time recorded by the behavior control system
        
        % Match timestamps -> in case of mismatch, try to fix
        if ~ismatch(ts,son2)
            son2 = clearttls(son2); % eliminate recorded TTL's within 0.5s from each other - broken TTL pulse
            if ~ismatch(ts,son2)
                son2 = trytomatch(ts,son2);  % try to match time series by shifting
                if ~ismatch(ts,son2)
                    [ts, keepinx] = trytomatch2(ts,son2);  % try to match time series by shifting
                    TE2 = shortenTE(TE2,keepinx);
                    if ~ismatch(ts,son2)
                        son2 = tryinterp(ts,son2); % interpolate missing TTL's or delete superfluous TTL's up to 10 erroneous TTl's
                        if ~ismatch(ts,son2)  % TTL matching failure
                            error('MakeTrialEvents:TTLmatch','Matching TTLs failed.')
                        else
                            warning('MakeTrialEvents:TTLmatch','Missing TTL interpolated.')
                        end
                    else
                        warning('MakeTrialEvents:TTLmatch','Shifted TTL series.')
                    end
                else
                    warning('MakeTrialEvents:TTLmatch','Shifted TTL series.')
                end
            else
                warning('MakeTrialEvents:TTLmatch','Broken TTLs cleared.')
            end
        end
        
        % Eliminate last TTL's recorded in only one system
        sto = TE2.StimulusOn;
        if length(son2) > length(ts)   % time not saved in behavior file (likely reason: autosave was used)
            son2 = son2(1:length(ts));
        elseif length(son2) < length(ts)  % time not recorded on Neuralynx (likely reason: recording stopped)
            shinx = 1:length(son2);
            ts    = ts(shinx);
            sto   = sto(shinx);
            TE2   = shortenTE(TE2,shinx);
        end
        TE2.TrialStart = son2 - sto;
        
        % Save synchronized 'TE_recording' file
        if ~isempty(TE2.TrialStart)
            save([sessionpath '\' 'TE_recording.mat'],'-struct','TE2')
        else
            error('MakeTrialEvents:noOutput','Synchronization process failed.');
        end
        
    case 'stim'
        BurstSeparation = 0.5;
        
        % Load Open Ephys events file -------------------------------------
        try
            load([sessionpath '\' 'TTLOnOff.mat'],'pulseon','pulseoff');
            A = dir([sessionpath '\' 'BA*']);
            A = A.name;
            [~, timestamps, ~] = load_open_ephys_data([sessionpath '\' A...
                '\all_channels.events']);
        catch
            A = dir([sessionpath '\' 'BA*']);
            A = A.name;
            [~,pulseon, pulseoff] = load_events([sessionpath '\' A]);
            [~, timestamps, ~] = load_open_ephys_data([sessionpath '\' A...
                '\all_channels.events']);
        end
        % There are 2 lasers: CH1 = Blue; CH2 = Yellow
        pulseID   = [ones(1,length(pulseon.CH1)), repmat(2,1,length(pulseon.CH2))]; 
        pulseon   = [pulseon.CH1; pulseon.CH2]'; 
        pulseoff  = [pulseoff.CH1; pulseoff.CH2]';
        nvalid_pulses = length(pulseon);  % calculating number of pulses
                
        % Preallocate stimulus events -------------------------------------
        SE = struct;
        SE.ProtocolName  = cell(1,nvalid_pulses); % e.g. LaserStimProtocol2
        SE.ProtocolID    = cell(1,nvalid_pulses); % 'S' for stim. and 'B' for behav. protocols
        SE.ProtocolStart = nan(1,nvalid_pulses);  % time stamp for protocol begin
        SE.ProtocolEnd   = nan(1,nvalid_pulses);  % time stamp for protocol end
        SE.StimType      = nan(1,nvalid_pulses);  % 1 for single pulse, 2 for burst stimualtion
        SE.TrialStart    = zeros(1,nvalid_pulses);% should be 0
        
        SE.BurstOn  = nan(1,nvalid_pulses);      % onset of first pulse in the burst %csak az elejenel
        SE.BurstOff = nan(1,nvalid_pulses);      % offset of last pulse in the burst
        SE.BurstDur = nan(1,nvalid_pulses);      % duration of burst (between last pulse offset and first pulse onset
        SE.BurstIBI = nan(1,nvalid_pulses);      % inter-burst interval
        SE.BurstNPulse  = nan(1,nvalid_pulses);  % number of pulses in burst (Events_EventStrings)
        SE.PrevBurstOff = nan(1,nvalid_pulses);  % end of previous burst
        SE.NextBurstOn  = nan(1,nvalid_pulses);  % start of next burst
        SE.PreBurstIBI  = nan(1,nvalid_pulses);  % time since PrevBurstOff
        SE.PostBurstIBI = nan(1,nvalid_pulses);  % time to NextBurstOn
%         SE.BurstID  = nan(1,nvalid_pulses);      % burst rank in the protocol (NTrial)
        
        SE.PulseID  = nan(1,nvalid_pulses);      % laser used: CH1 = Blue; CH2 = Yellow
        SE.PulseOn  = nan(1,nvalid_pulses);      % onset of stim. pulse
        SE.PulseOff = nan(1,nvalid_pulses);      % offset of stim. pulse
        SE.PulseDur = nan(1,nvalid_pulses);      % duration of light pulse (Events_EventStrings)
        SE.PulseIPI = nan(1,nvalid_pulses);      % inter-pulse interval (Events_EventStrings)
        SE.PulseFreq    = nan(1,nvalid_pulses);  % pulse frequency (Events_EventStrings)
%         SE.PulsePower   = nan(1,nvalid_pulses);  % stimulus intensity (Events_EventStrings)
        SE.PrevPulseOff = nan(1,nvalid_pulses);  % offset of previous pulse
        SE.NextPulseOn  = nan(1,nvalid_pulses);  % onset of next pulse
        SE.PrePulseIPI  = nan(1,nvalid_pulses);  % time since PrevPulseOff
        SE.PostPulseIPI = nan(1,nvalid_pulses);  % time to NextPulseOn
        
%         SE.FirstPulse = nan(1,nvalid_pulses);    % is this the first pulse in a burst?
        SE.ZeroPulse  = nan(1,nvalid_pulses);    % extrapolate one pulse time back
%         SE.LastPulse  = nan(1,nvalid_pulses);    % is this the last pulse in a burst?
%         SE.PulseNum   = nan(1,nvalid_pulses);    % pulse rank in the burst
        
        % Filling SE struct -----------------------------------------------
        SE.PulseID  = pulseID;
        SE.PulseOn  = pulseon;
        SE.PulseOff = pulseoff;
        SE.ProtocolStart = repmat(timestamps(3),1,nvalid_pulses);
        SE.ProtocolEnd   = repmat(timestamps(end),1,nvalid_pulses); 
        SE.StimType      = repmat(2,1,nvalid_pulses); % 1 for single pulse, 2 for burst stimualtion
           
        for currentEvent = 1:nvalid_pulses
            SE.PulseDur(currentEvent)      = SE.PulseOff(currentEvent)-SE.PulseOn(currentEvent); % duration of stim. pulse (Events_EventStrings)
            SE.ProtocolName(currentEvent)  = cellstr('Tagging'); % e.g. LaserStimProtocol2
            SE.ProtocolID(currentEvent)    = cellstr('S');       % 'S' for stim. and 'B' for behav. protocol
        end
        
        for currentEvent = 2:nvalid_pulses
            SE.PulseIPI(currentEvent) = SE.PulseOn(currentEvent)-SE.PulseOff(currentEvent-1); % inter-pulse interval (Events_EventStrings)
        end
        
        % BurstOn
        burstinx             = SE.PulseIPI > BurstSeparation;
        SE.BurstOn(burstinx) = SE.PulseOn(burstinx);
        SE.BurstOn(1)        = SE.PulseOn(1);
        
        % BurstOff
        burstinx2              = find(burstinx) - 1;
        SE.BurstOff(burstinx2) = SE.PulseOff(burstinx2);
        SE.BurstOff(end)       = SE.PulseOff(end);
        
        % Inter-burst interval
        SE.PrevBurstOff = [SE.ProtocolStart(1) SE.BurstOff(1:end-1)];   % previous burst offset (first protocol start for first pulse)
        SE.NextBurstOn  = [SE.BurstOn(2:end) SE.ProtocolEnd(end)];   % next burst onset (last protocol end for last pulse)
        SE.BurstIBI     = SE.NextBurstOn - SE.BurstOff; % inter-burst interval
        SE.PreBurstIBI  = SE.BurstOn - SE.PrevBurstOff; % previous inter-burst interval
        SE.PostBurstIBI = SE.NextBurstOn - SE.BurstOff; % next inter-burst interval
        
        % Inter-pulse interval
        SE.PrevPulseOff = [SE.ProtocolStart(1) SE.PulseOff(1:end-1)];   % previous pulse offset (first protocol start for first pulse)
        SE.NextPulseOn  = [SE.PulseOn(2:end) SE.ProtocolEnd(end)];   % next pulse onset (last protocol end for last pulse)
        SE.PrePulseIPI  = SE.PulseOn - SE.PrevPulseOff;   % inter-pulse interval before pulse
        SE.PostPulseIPI = SE.NextPulseOn - SE.PulseOff;   % inter-pulse interval after pulse
        
        % Burst parameters
        burston  = find(~isnan(SE.BurstOn));
        burstoff = find(~isnan(SE.BurstOff));
        for currentEvent = 1:nvalid_pulses
            prevburston  = burston(find(burston<=currentEvent,1,'last'));
            nextburstoff = burstoff(find(burstoff>=currentEvent,1,'first'));
            SE.BurstNPulse(currentEvent) = nextburstoff - prevburston + 1; % number of pulses in burst
            SE.BurstDur(currentEvent)    = SE.PulseOn(nextburstoff) - SE.PulseOn(prevburston);
            SE.PulseFreq(currentEvent)   = SE.BurstNPulse(currentEvent) / SE.BurstDur(currentEvent);
        end
        
        % Save
        save([sessionpath '\' 'TE_StimEvent.mat'], '-struct','SE'); %Save OE timestamps into a .mat file
end

% -------------------------------------------------------------------------
function I = ismatch(ts,son2)

% Check if the two time series match notwithstanding a constant drift
clen = min(length(ts),length(son2));
I = abs(max(diff(ts(1:clen)-son2(1:clen)))) < 0.05;  % the difference between the timestamps on 2 systems may have a constant drift, but it's derivative should still be ~0

% note: abs o max is OK, the derivative is usually a small neg. number due
% to drift of the timestamps; max o abs would require a higher tolerance
% taking the drift into account (if 2 event time stamps are far, the drift
% between them can be large)

% -------------------------------------------------------------------------
function son2 = tryinterp(ts,son2)

% Interpolate missing TTL's or delete superfluous TTL's up to 10 erroneous TTl's
for k = 1:10
    if ~ismatch(ts,son2)
        son3 = son2 - son2(1) + ts(1);
        adt = diff(ts(1:min(length(ts),length(son2)))-son2(1:min(length(ts),length(son2))));
        badinx = find(abs(adt)>0.1,1,'first') + 1;  % find problematic index
        if adt(badinx-1) < 0    % interploate
            ins = ts(badinx) - linterp([ts(badinx-1) ts(badinx+1)],[ts(badinx-1)-son3(badinx-1) ts(badinx+1)-son3(badinx)],ts(badinx));
            son2 = [son2(1:badinx-1) ins+son2(1)-ts(1) son2(badinx:end)];
        else
%             ins = son3(badinx) - linterp([son3(badinx-1) son3(badinx+1)],[son3(badinx-1)-ts(badinx-1) son3(badinx+1)-ts(badinx)],son3(badinx));
%             ts = [ts(1:badinx-1) ins ts(badinx:end)];
            son2(badinx) = [];   % delete
        end
    end
end

% -------------------------------------------------------------------------
function son2 = trytomatch(ts,son2)

% Try to match time series by shifting
len = length(son2) - 15;
minx = nan(1,len);
for k = 1:len
    minx(k) = max(diff(ts(1:15)-son2(k:k+14)));  % calculate difference in the function of shift
end
mn = min(abs(minx));
minx2 = find(abs(minx)==mn);
minx2 = minx2(1);   % find minimal difference = optimal shift
if mn < 0.2
    son2_temp = son2(minx2:min(minx2+length(ts)-1,length(son2)));
    if length (son2_temp) < (len-10)
        warning ('Lost too many TTL-s.')
    else
        son2 = son2_temp;
    end
    
end

% -------------------------------------------------------------------------
function [ts, keepinx] = trytomatch2(ts,son2)

% Try to match time series by shifting
keepinx = [];
len = length(son2) - 15;
minx = nan(1,len);
for k = 1:len
    minx(k) = max(diff(son2(1:15)-ts(k:k+14)));  % calculate difference in the function of shift
end
mn = min(abs(minx));
minx2 = find(abs(minx)==mn);
minx2 = minx2(1);   % find minimal difference = optimal shift
if mn < 0.2
    keepinx = minx2:length(ts);
    %     ts = ts(minx2:min(minx2+length(son2)-1,length(ts)));
    ts_temp = ts(keepinx);
    if length (ts_temp) < (len-10)
        warning ('Lost too many TTL-s.')
    else
        ts = ts_temp;
    end
    
end

% -------------------------------------------------------------------------
function son2 = clearttls(son2)

% Eliminate recorded TTL's within 0.5s from each other
inx = [];
for k = 1:length(son2)-1
    s1 = son2(k);
    s2 = son2(k+1);
    if s2 - s1 < 0.5
        inx = [inx k+1]; %#ok<AGROW>
    end
end
son2(inx) = [];

% -------------------------------------------------------------------------
function TE2 = shortenTE(TE2,shinx)

% Eliminate behavioral trials
fnm = fieldnames(TE2);
for k = 1:length(fieldnames(TE2))
    TE2.(fnm{k}) = TE2.(fnm{k})(shinx);
end
