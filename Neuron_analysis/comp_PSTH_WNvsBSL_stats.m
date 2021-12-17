
function stats = comp_PSTH_WNvsBSL_stats(cellid,event,baselinewin,testwin,varargin) %#ok<*INUSL>
%comp_PSTH_WNvsBSL_stats   Peri-stimulus time histogram stats.
%  Given PSTH activity in baseline and a test window, the functions calculates
%  if they are significantly different (with plots) 
%
%   Mandatory input arguments:
%   CELLID: defines the cell 
%   EVENT: the event to which the PSTH windows are aligned: the first
%    event is used for the baseline, the second for the test
%   BASELINEWIN: window for calculation relative to the baseline window in seconds
%   TESTWIN: window for calculation relative to the test window in seconds
%
% -------------------------------------------------------------------------
% Based on ULTIMATE_PSTH_WM_VS_BSLN by Nicola Solari
% 
% Cecília Pardo-Bellver, 2021
% Laboratory of Network Neurophysiology
% Institute of Experimental Medicine, Hungary.
% ------------------------------------------------------------------------- 

% Default arguments
prs = inputParser;
addRequired(prs,'cellid')
addRequired(prs,'event',@(s)isnumeric(s)|(iscellstr(s)&isequal(length(s),2)))   % reference event
addRequired(prs,'baselinewin',@(s)isnumeric(s)&isequal(length(s),2))  % time window relative to the event for stat. testing, in seconds
addRequired(prs,'testwin',@(s)isnumeric(s)&isequal(length(s),2))  % time window relative to the event for stat. testing, in seconds

addParameter(prs,'cellpath','L:\_Neurons');
addParameter(prs,'filterinput',[])   % some filters need additional input
addParameter(prs,'maxtrialno',5000)   % downsample events if more than 'maxtrialno'
addParameter(prs,'first_event','',@(s)isempty(s)|ischar(s))   % exclude spikes before previous event
addParameter(prs,'last_event','',@(s)isempty(s)|ischar(s))   % exclude spikes after following events
addParameter(prs,'dt',0.001,@isnumeric)   % time resolution of the binraster, in seconds
addParameter(prs,'sigma',0.02,@isnumeric)     % smoothing kernel for the smoothed PSTH
addParameter(prs,'margin',[-0.1 0.1])  % margins for PSTH calculation to get rid of edge effect due to smoothing
addParameter(prs,'relative_threshold',0.5,@(s)isnumeric(s)&s>=-1&s<=1)   % threshold used to assess interval limits in PSTH_STATS; negative thresholds selects the full window
addParameter(prs,'sig_thr',0.01) % significance threshold 
addParameter(prs,'window',[-3 3]) % significance threshold 
addParameter(prs,'tag',{}) 
parse(prs,cellid,event,baselinewin,testwin,varargin{:})
g = prs.Results;

% Load spikes
load([cellid '.mat'],'TS')

times = g.window(1):g.dt:g.window(2);
times = round(times,4);
stimes       = ((TS(:,1)/10000));
event_spikes = cell(1,length(g.event));
spk    = zeros(length(g.event), length(times));
for nTT = 1:length(g.event)
    event_spikes{nTT} = stimes(stimes>=g.event(nTT)+g.window(1) & ...
        stimes<g.event(nTT)+g.window(2)) - g.event(nTT);
    ind_spikes        = round((event_spikes{nTT}-times(1))/g.dt) + 1;
    if ~isempty(ind_spikes)
        spk(nTT,ind_spikes) = 1; 
    end
end
        
% Select spikes, calculate bin raster
time_b = round((g.baselinewin(1)+g.margin(1)):g.dt:(g.baselinewin(2)+g.margin(2)),4);
time_t = round((g.testwin(1)+g.margin(1)):g.dt:(g.testwin(2)+g.margin(2)),4);

spk_b = spk(:,ismember(times,time_b)==1);
spk_t = spk(:,ismember(times,time_t)==1);

% PSTH & SPSTH ------------------------------------------------------------
psth_b  = mean(spk_b)/g.dt;
psth_t  = mean(spk_t)/g.dt;

kernel = gausswin(2/g.sigma); % Create a gaussian as kernel for the convolution
conv_psth_b = conv(psth_b,kernel); % Convolution of the psth
last_b  = length(conv_psth_b) - floor(length(kernel)/2); % Cutting excess time
first_b = 1 + last_b - length(psth_b);
spsth_b = conv_psth_b(first_b:last_b)/sum(kernel);
conv_psth_t = conv(psth_t,kernel);
last_t  = length(conv_psth_t) - floor(length(kernel)/2);
first_t = 1 + last_t - length(psth_t);
spsth_t = conv_psth_t(first_t:last_t)/sum(kernel);

% Readjust points of the psth
time0_b = abs(g.baselinewin(1)+g.margin(1)) * (1 / g.dt); % zero point (diveding directly with 'g.dt' can result in numeric deviation from the desired integer result)
time_b  = round(time0_b); % still numeric issues
if abs(time_b-time0_b) > 1e-10
    error('Zero point is not an integer.')
end
inx_b   = (time_b+1+g.baselinewin(1)/g.dt):(time_b+1+g.baselinewin(2)/g.dt); % indices for cutting margin
psth_b  = psth_b(:,inx_b);
spsth_b = spsth_b(:,inx_b);

time0_t = abs(g.testwin(1)+g.margin(1)) * (1 / g.dt);   % zero point (diveding directly with 'g.dt' can result in numeric deviation from the desired integer result)
time_t  = round(time0_t);   % still numeric issues
if abs(time_t-time0_t) > 1e-10
    error('Zero point is not an integer.')
end
inx_t   = (time_t+1+g.testwin(1)/g.dt):(time_t+1+g.testwin(2)/g.dt); % indices for cutting margin
psth_t  = psth_t(:,inx_t);
spsth_t = spsth_t(:,inx_t);

% Output statistics -------------------------------------------------------
statpsth_b = spsth_b;   % use smoothed PSTH for finding activation and inhibition windows in psth_stats
statpsth_t = spsth_t;

figure
stats = comp_PSTH_WNvsBSL_statszero(1,spk_b,spk_t,...
    statpsth_b,statpsth_t,g.dt,g.baselinewin, g.testwin, g.tag, g.sig_thr);
