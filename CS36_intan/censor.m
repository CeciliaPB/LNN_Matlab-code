function allspikes = censor(allspikes,sr,deadtime)
%CENSOR   Apply censoring period after detected peaks.
%   ALLSPIKES = CENSOR(ALLSPIKES,SR,DEADTIME) deletes points of the
%   ALLSPIKES series sampled at SR sampling rate that are within a
%   censoring period from a previous point. The censoring period is defined
%   by DEADTIME, given in seconds.
%
%   Note that as CENSOR was meant to be a low level utility function,
%   ALLSPIKES should be given in data points (not converted into seconds).
%
%   See also PEAKDETECT and CUTSPIKE.

%   Balazs Hangya, TENSS 2016
%   hangya.balazs@koki.mta.hu

% Recording dead time
dtp = deadtime * sr;   % dead time in data points

% Apply censoring
isi = diff(allspikes);   % inter-spike intervals
censor_inx = isi < dtp;   % ISI < dead time
allspikes([false censor_inx]) = [];   % remove spikes within the censor period