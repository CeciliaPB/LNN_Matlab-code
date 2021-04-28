function [events,epochs] = defineEventsEpochs_gambling
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

% Define events and epochs
%              EventLabel           EventTrigger1        EventTrigger2        Window
i = 1;
events(i,:) = {'StimulusOn',        'StimulusOn',        'StimulusOn',        [-6 6]};    i = i + 1;
events(i,:) = {'StimulusOff',       'StimulusOff',       'StimulusOff',       [-6 6]};    i = i + 1;
events(i,:) = {'DeliverFeedback',   'DeliverFeedback',   'DeliverFeedback',   [-6 6]};    i = i + 1;
events(i,:) = {'DeliverAllFeedback','DeliverAllFeedback','DeliverAllFeedback',[-6 6]};    i = i + 1;
events(i,:) = {'LickLIn',           'LickLIn',            'LickLIn',          [-6 6]};    i = i + 1;
events(i,:) = {'LickLOut',          'LickLOut',           'LickLOut',         [-6 6]};    i = i + 1;
events(i,:) = {'LickRIn',           'LickRIn',            'LickRIn',          [-6 6]};    i = i + 1;
events(i,:) = {'LickROut',          'LickROut',           'LickROut',         [-6 6]};    i = i + 1;

% Variable events
% events(i,:) = {'StimulusSampling','StimulusOn',    'StimulusOff',     [-6 6]};    i = i + 1;

% Define epochs for rate calculations
%               EpochLabel      ReferenceEvent      FixedWindow       RealWindow
i = 1;
epochs(i,:) = {'StimulusOn',    'StimulusOn',       [0.0 0.3],        'StimulusSampling'};    i = i + 1;