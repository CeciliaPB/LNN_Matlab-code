function [events,epochs] = defineEventsEpochs_gambling(varargin)
%DEFINEEVENTSEPOCHS_GAMBLING   Define events and epochs for spike extraction.
%   [EVENTS,EPOCHS] = DEFINEEVENTSEPOCHS_GAMBLING defines events and epochs
%   for spike extraction. 
%
%   EVENTS is a Nx4 cell array with columns corresponding to EventLabel,
%   EventTrigger1, EventTrigger2, Window. EventLabel is the name for
%   referencing the event. EventTrigger1 and EventTrigger2 are names of
%   TrialEvent variables (e.g. 'LeftPortIn'). For fixed windows, the two
%   events are the same; for variable windows, they correspond to the start
%   and end events. Window specifies time offsets relative to the events;
%   e.g. events(1,:) = {'OdorValveOn','OdorValveOn','OdorValveOn',[-3 3]};
%
%   EPOCH is a Nx4 cell array with columns corresponding to  EpochLabel, 
%   ReferenceEvent, Window, RealWindow. EventLabel is the name for 
%   referencing the epoch. ReferenceEvent should match an EventLabel in 
%   EVENTS (used for calculating the epoch rates). RealWindow is currently
%   not implemented (allocated for later versions).
%
%   DEFINEEVENTSEPOCHS_GAMBLING defines events and epochs for auditory
%   two lickports gambling task.
%
%   See also MAKETRIALEVENTS2_GONOGO and DEFINEEVENTSEPOCHS_DEFAULT.
% -------------------------------------------------------------------------
%   Cecília Pardo-Bellver, 2021
%   Laboratory of Network Neurophysiology
%   Institute of Experimental Medicine, Hungary.
%
%   Based on DEFINEEVENTSEPOCHS_DEFAULT by Balazs Hangya. Edit log: BH 5/2/14
% -------------------------------------------------------------------------

eventtype = 'behav';
if nargin == 1
    eventtype = varargin{1};
end

switch  eventtype
    case {'behav','event'}
        % Define events and epochs
        %              EventLabel           EventTrigger1        EventTrigger2        Window
        ii = 1;
        events(ii,:) = {'StimulusOn',        'StimulusOn',        'StimulusOn',        [-6 6]};    ii = ii + 1;
        events(ii,:) = {'StimulusOff',       'StimulusOff',       'StimulusOff',       [-6 6]};    ii = ii + 1;
        events(ii,:) = {'DeliverFeedback',   'DeliverFeedback',   'DeliverFeedback',   [-6 6]};    ii = ii + 1;
        events(ii,:) = {'DeliverAllFeedback','DeliverAllFeedback','DeliverAllFeedback',[-6 6]};    ii = ii + 1;
        events(ii,:) = {'LickLIn',           'LickLIn',            'LickLIn',          [-6 6]};    ii = ii + 1;
        events(ii,:) = {'LickLOut',          'LickLOut',           'LickLOut',         [-6 6]};    ii = ii + 1;
        events(ii,:) = {'LickRIn',           'LickRIn',            'LickRIn',          [-6 6]};    ii = ii + 1;
        events(ii,:) = {'LickROut',          'LickROut',           'LickROut',         [-6 6]};

        % Variable events
        % events(i,:) = {'StimulusSampling','StimulusOn',    'StimulusOff',     [-6 6]};    i = i + 1;

        % Define epochs for rate calculations
        %               EpochLabel      ReferenceEvent      FixedWindow       RealWindow
        kk = 1;
        epochs(kk,:) = {'StimulusOn',    'StimulusOn',       [0.0 0.3],        'StimulusSampling'};
    
    case 'stim'
        %              EventLabel       EventTrigger1    EventTrigger2  Window
        ii = 1;
        events(ii,:) = {'BurstOn',       'BurstOn',      'BurstOn',      [-6 6]};   ii = ii + 1;
        events(ii,:) = {'PulseOn',       'PulseOn',      'PulseOn',      [-3 3]};   ii = ii + 1;
        events(ii,:) = {'PreBurstIBI',   'PrevBurstOff', 'BurstOn',      [0 0]};    ii = ii + 1;
        events(ii,:) = {'BurstPeriod',   'BurstOn',      'BurstOff',     [0 0]};    ii = ii + 1;
        events(ii,:) = {'NextBurstIBI',  'BurstOff',     'NextBurstOn',  [0 0]};    ii = ii + 1;

        % Variable events
        events(ii,:) = {'PreBurstIBI2',  'PrevBurstOff', 'BurstOn',      [-1 0]};   ii = ii + 1;
        events(ii,:) = {'BurstPeriod2',  'BurstOn',      'BurstOff',     [-1 0]};   ii = ii + 1;
        events(ii,:) = {'NextBurstIBI2', 'BurstOff',     'NextBurstOn',  [-1 0]};   ii = ii + 1;
        % events(i,:) = {'OmitPulse',     'OmitPulse',    'OmitPulse',    [-6 6]};   i = i + 1;
        events(ii,:) = {'ZeroPulse',     'ZeroPulse',    'ZeroPulse',    [-6 6]};

        % Define epochs for rate calculations
        %               EpochLabel             ReferenceEvent  FixedWindow          RealWindow
        kk = 1;
        epochs(kk,:) = {'FixedBaseline5',       'BurstOn',      [-0.005 0.0],  'PreBurstIBI'};   kk = kk + 1;
        epochs(kk,:) = {'FixedLightResponse5',  'BurstOn',      [0.0 0.005],   'BurstPeriod'};   kk = kk + 1;
        epochs(kk,:) = {'FixedBaseline5a',      'BurstOn',      [-0.005 0.0],  'BurstOn'};       kk = kk + 1;
        epochs(kk,:) = {'FixedLightResponse5a', 'BurstOn',      [0.0 0.005],   'BurstOn'};       kk = kk + 1;
        epochs(kk,:) = {'FixedBaseline10',      'BurstOn',      [-0.01 0.0],   'PreBurstIBI'};   kk = kk + 1;
        epochs(kk,:) = {'FixedLightResponse10', 'BurstOn',      [0.0 0.01],    'BurstPeriod'};   kk = kk + 1;
        epochs(kk,:) = {'FixedBaseline20',      'BurstOn',      [-0.02 0.0],   'PreBurstIBI'};   kk = kk + 1;
        epochs(kk,:) = {'FixedLightResponse20', 'BurstOn',      [0.0 0.02],    'BurstPeriod'};   kk = kk + 1;
        epochs(kk,:) = {'FixedBaseline40',      'BurstOn',      [-0.04 0.0],   'PreBurstIBI'};   kk = kk + 1;
        epochs(kk,:) = {'FixedLightResponse40', 'BurstOn',      [0.0 0.04],    'BurstPeriod'};   kk = kk + 1;
        epochs(kk,:) = {'BurstBaseline',        'PreBurstIBI',  [NaN NaN],     'NaN'};           kk = kk + 1;
        epochs(kk,:) = {'BurstResponse',        'BurstPeriod',  [NaN NaN],     'NaN'};

end

end