function [iFR,tFR] = instFiringRate(timestamps,time,bins)

% instFiringRate - Instantaneous firing rate of a spike train
%   iFR = instFiringRate(timestamps, time) calculates the instantaneous
%   firing rate at times: tFR of a spike train with spikes at times:
%   timestamps. In the interval between two spikes n and n+1, the
%   instanteous firing rate: iFR is defined as:
%
%      iFR = 1 / (timestamps(n+1) - timestamps(n).
%
%   Even though there is no natural definition for the iFR before the first
%   spike or after the last spike in a train, this function sets iFR to zero
%   in those intervals.
%
% INPUTS: 
%   - timestamps: vector with timestamps of spike events.
%   - time: vector, only needs start and end times. Make sure the units in
%   time and timestamps are matching. 
%   - bins:  bin size of interval. If not provided: Default 100 = 0.01s.
%
% OUTPUTS: 
%   - iFR: Instantaneous firing rate. Vector.
%   - tFR: Calculated time vector.
%
% Examples 
%   [iFR,tFR] = instFiringRate(timestamps, time); 
%   [iFR,~]   = instFiringRate(timestamps,time,100);
%
% -------------------------------------------------------------------------
% Based on instantfr.m by Daniel A. Wagenaar.
% 
% Cecília Pardo-Bellver, 2020
% Laboratory of Network Neurophysiology
% Instinute of Experimental Medicine, Hungary.
% -------------------------------------------------------------------------
 
% Params ------------------------------------------------------------------
epsi = 1e-6;
if nargin < 3
    bins = 100; % Default sampling: 0.01 s.
    bins = 1/bins; 
    tFR = min(time):bins:max(time);
else
    bins = 1/bins; 
    tFR = min(time):bins:max(time);
end

% Calculate time difference between spikes --------------------------------
tt0 = [timestamps(:)' - epsi; timestamps(:)' + epsi];
idt = [diff([0 timestamps(:)']); diff([timestamps(:)' inf])];
tt0 = tt0(:)';
idt = idt(:)';

% Calculate the instantaneous firing rate ---------------------------------
if nargin < 2
  iFR = 1./idt(2:end-1);
  tFR = tt0(2:end-1);
  if size(timestamps,1)~=1
    iFR=iFR';
    tFR=tFR';
  end
else
  if isempty(timestamps)
    iFR=0*tFR;
    return;
  end
  idt1 = inf;%tt0(1)-tt(1);
  idtn = inf;%tt(end)-tt0(end);
  ifr0 = 1./[idt1 idt1 idt(2:end-1) idtn idtn];
  tt0 = [tFR(1) tt0 tFR(end)];
  iFR = interp1(tt0,ifr0,tFR,'linear');
  if size(tFR,1)~=1
    iFR=iFR';
  end
end

if nargout<2
  clear tt
end

end