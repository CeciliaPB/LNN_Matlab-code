function [spt, pk] = disc(unit,thres)
%DISC   Unit discrimination.
%   [SPT PK] = DISC(UNIT,THR) returns the location (SPT) and value (PK) of
%   UNIT peaks above THR.

% Input argument check
narginchk(2,2)
unit = unit(:)';   % row vector

% Segments above threshold
dsc = find(unit>=thres);   % find all points above threshold 
if isempty(dsc)  % no spikes
    spt = [];
    pk = [];
    return
end
lendsc = length(dsc);
df = diff(dsc);   % distance of points above threshold
gaps = find(df>1);   % find gaps between segments above threshold
gaps = [0 gaps lendsc];   % there's a spike between every two gaps

% Find peaks for each segment
NumSpk = length(gaps) - 1;   % number of spikes
[spt, pk] = deal(nan(1,NumSpk));
for iS = 1:NumSpk
    [mx, mxloc] = max(unit(dsc(gaps(iS)+1:gaps(iS+1))));
    spt(iS) = dsc(1,gaps(iS)+mxloc);
    pk(iS) = mx;   % value at peak
end