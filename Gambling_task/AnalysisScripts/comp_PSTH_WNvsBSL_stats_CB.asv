
function stats = comp_PSTH_WNvsBSL_stats_CB(cellid,event,baselinewin,testwin,varargin) %#ok<*INUSL>
%comp_PSTH_WNvsBSL_stats_CB   Peri-stimulus time histogram stats.
%  Given PSTH activity in baseline and a test window, the functions calculates
%  if they are significantly different (with plots) 
%
%   Mandatory input arguments:
%   CELLID: defines the cell (see CellBase documentation) 
%   EVENT: the event to which the PSTH windows are aligned: the first
%    event is used for the baseline, the second for the test
%   BASELINEWIN: window for calculation relative to the baseline window in seconds
%   TESTWIN: window for calculation relative to the test window in seconds
%
%   See also comp_PSTH_WNvsBSL_main, ULTIMATE_PSTH, PSTH_STATS, STIMES2BINRASTER, BINRASTER2PSTH, BINRASTER2APSTH,
%   APSTH, VIEWCELL2B, PARTITION_TRIALS and FILTERTRIALS.
% -------------------------------------------------------------------------
% Based on ULTIMATE_PSTH_WM_VS_BSLN by Nicola Solari
% 
% Cecília Pardo-Bellver, 2021
% Laboratory of Network Neurophysiology
% Institute of Experimental Medicine, Hungary.
% ------------------------------------------------------------------------- 

% Default arguments
prs = inputParser;
addRequired(prs,'cellid',@(s)iscellid(s)|issessionid(s))
addRequired(prs,'event',@(s)ischar(s)|(iscellstr(s)&isequal(length(s),2))|...
    isa(s,'function_handle'))   % reference event
addRequired(prs,'baselinewin',@(s)isnumeric(s)&isequal(length(s),2))  % time window relative to the event for stat. testing, in seconds
addRequired(prs,'testwin',@(s)isnumeric(s)&isequal(length(s),2))  % time window relative to the event for stat. testing, in seconds

addParameter(prs,'cellpath','L:\Cecilia\Behaviour\Risk assesment');
addParameter(prs,'event_filter','none',@(s)ischar(s)|iscellstr(s))   % filter events based on properties
addParameter(prs,'filterinput',[])   % some filters need additional input
addParameter(prs,'maxtrialno',5000)   % downsample events if more than 'maxtrialno'
addParameter(prs,'first_event','',@(s)isempty(s)|ischar(s))   % exclude spikes before previous event
addParameter(prs,'last_event','',@(s)isempty(s)|ischar(s))   % exclude spikes after following events
addParameter(prs,'dt',0.001,@isnumeric)   % time resolution of the binraster, in seconds
addParameter(prs,'sigma',0.02,@isnumeric)     % smoothing kernel for the smoothed PSTH
addParameter(prs,'margin',[-0.1 0.1])  % margins for PSTH calculation to get rid of edge effect due to smoothing
addParameter(prs,'partitions','all')   % partition trials
addParameter(prs,'isadaptive',true,@(s)islogical(s)|ismember(s,[0 1 2]))   % use adaptive PSTH algorithm
addParameter(prs,'relative_threshold',0.5,@(s)isnumeric(s)&s>=-1&s<=1)   % threshold used to assess interval limits in PSTH_STATS; negative thresholds selects the full window
addParameter(prs,'sig_thr',0.01) % significance threshold 
addParameter(prs,'window',[-3 3]) % significance threshold 
addParameter(prs,'display',false,@(s)islogical(s)|ismember(s,[0 1]))   % control displaying rasters and PSTHs
parse(prs,cellid,event,baselinewin,testwin,varargin{:})
g = prs.Results;

% Load trial events
try
    TE = load([g.cellpath '\' cellid(1:5) '\' cellid(7:15) '\' 'TE_recording.mat']);   % load events
    SP = loadcb(cellid,'EVENTSPIKES');   % load prealigned spikes
catch ME
    disp('There was no behavioral protocol for this session.')
    error(ME.message)
end
        
% Events, time, valid trials
time_b = (g.baselinewin(1)+g.margin(1)):g.dt:(g.baselinewin(2)+g.margin(2));
time_t = (g.testwin(1)+g.margin(1)):g.dt:(g.testwin(2)+g.margin(2));

% Different ref. event for test and baseline window
event1 = event{1};   % baseline event
event2 = event{2};   % test event
event_pos1 = findcellstr(SP.events(:,1),event1);
event_pos2 = findcellstr(SP.events(:,1),event2);
if event_pos1 * event_pos2 == 0
    error('Event name not found');
end
stimes1 = SP.event_stimes{event_pos1};   % baseline spike times
stimes2 = SP.event_stimes{event_pos2};   % test spike times
triggerevent1 = SP.events{event_pos1,2};   % trigger event for baseline (event name may differ)
triggerevent2 = SP.events{event_pos2,2};   % trigger event for test (event name may differ)

valid_trials1 = 1:length(stimes1);
valid_trials2 = 1:length(stimes1);
if ~isequal(valid_trials1,valid_trials2)
    error('Valid trials should be the same for baseline and test period.')
else
    valid_trials = valid_trials1;
end

if ~isempty(g.first_event) || ~isempty(g.last_event)
    [stimes1, starttimes1, endtimes1] = dynamicSpikeWindow(stimes1,TE,...
        triggerevent1,g.first_event,g.last_event); % restrict between previous and following event
    [stimes2, starttimes2, endtimes2] = dynamicSpikeWindow(stimes2,TE,...
        triggerevent2,g.first_event,g.last_event); % restrict between previous and following event
end

% Calculate bin rasters
% Different ref. event for test and baseline window
spk_b = stimes2binraster(stimes1,time_b,g.dt);
spk_t = stimes2binraster(stimes2,time_t,g.dt);
if ~isempty(g.first_event) || ~isempty(g.last_event)
    spk_b = nanpadspk(time_b,spk_b,starttimes1,endtimes1);  % replace zeros with NaNs outside the dynamic window
    spk_t = nanpadspk(time_t,spk_t,starttimes2,endtimes2); 
end

% Partition trials
[COMPTRIALS, tags] = partition_trials(TE,g.partitions);

% PSTH --------------------------------------------------------------------
switch g.isadaptive
    case {0,false}
        [psth_b, spsth_b, ~] = binraster2psth(spk_b,g.dt,g.sigma,...
            COMPTRIALS,valid_trials);
        [psth_t, spsth_t, ~] = binraster2psth(spk_t,g.dt,g.sigma,...
            COMPTRIALS,valid_trials);
    case {1, true}
        [psth_b, spsth_b, ~] = binraster2apsth(spk_b,g.dt,g.sigma,...
            COMPTRIALS,valid_trials);
        [psth_t, spsth_t, ~] = binraster2apsth(spk_t,g.dt,g.sigma,...
            COMPTRIALS,valid_trials);        
    case 2
        [psth_b, spsth_b, ~] = binraster2dapsth(spk_b,g.dt,g.sigma,...
            COMPTRIALS,valid_trials);
        [psth_t, spsth_t, ~] = binraster2dapsth(spk_t,g.dt,g.sigma,...
            COMPTRIALS,valid_trials);        
end

% Readjust points of the psth
time0_b = abs(g.baselinewin(1)+g.margin(1)) * (1 / g.dt); % zero point (diveding directly with 'g.dt' can result in numeric deviation from the desired integer result)
time_b  = round(time0_b); % still numeric issues
if abs(time_b-time0_b) > 1e-10
    error('Zero point is not an integer.')
end
inx_b   = (time_b+1+g.baselinewin(1)/g.dt):(time_b+1+g.baselinewin(2)/g.dt); % indices for cutting margin
psth_b  = psth_b(:,inx_b);
spsth_b = spsth_b(:,inx_b);
NumPartitions = size(psth_b,1);
if NumPartitions > 1   % partitions
    pspk_b = spk_b;
    spk_b  = cell(1,NumPartitions);
    for iP = 1:NumPartitions
        spk_b{iP} = pspk_b(intersect(valid_trials,COMPTRIALS{iP}),inx_b);
    end
else
    spk_b = spk_b(valid_trials,inx_b);
end

time0_t = abs(g.testwin(1)+g.margin(1)) * (1 / g.dt);   % zero point (diveding directly with 'g.dt' can result in numeric deviation from the desired integer result)
time_t  = round(time0_t);   % still numeric issues
if abs(time_t-time0_t) > 1e-10
    error('Zero point is not an integer.')
end
inx_t   = (time_t+1+g.testwin(1)/g.dt):(time_t+1+g.testwin(2)/g.dt); % indices for cutting margin
psth_t  = psth_t(:,inx_t);
spsth_t = spsth_t(:,inx_t);
NumPartitions = size(psth_t,1);
if NumPartitions > 1   % partitions
    pspk_t = spk_t;
    spk_t  = cell(1,NumPartitions);
    for iP = 1:NumPartitions
        spk_t{iP} = pspk_t(intersect(valid_trials,COMPTRIALS{iP}),inx_t);
    end
else
    spk_t = spk_t(valid_trials,inx_t);
end

% Output statistics -------------------------------------------------------
switch g.isadaptive
    case {0,false}
        statpsth_b = spsth_b;   % use smoothed PSTH for finding activation and inhibition windows in psth_stats
        statpsth_t = spsth_t;
    case {1,2,true}
        statpsth_b = psth_t;   % use adaptive Spike Density Function for finding activation and inhibition windows in psth_stats
        statpsth_t = psth_t;
end

figure
if NumPartitions > 1
    stats   = cell(1,NumPartitions);   %#ok<PREALL> % return multiple stats if partitioning
    [stats] = comp_PSTH_WNvsBSL_statszero(NumPartitions,spk_b,spk_t,...
        statpsth_b,statpsth_t,g.dt,g.baselinewin, g.testwin, tags, g.sig_thr);
    for iP = 1:NumPartitions
        stats{iP}.tag = tags{iP};
    end
else
    %         this should not happen
end

% -------------------------------------------------------------------------
function spk = nanpadspk(time,spk,starttimes,endtimes)

% For variable windows, change padding to NaN to ensure correct averaging
NUMtrials = size(spk,1);   % number of trials
for iT = 1:NUMtrials    % loop through trials
    inx = time < starttimes(iT);
    spk(iT,inx) = NaN;   % NaN trials before previous event
end
for iT = 1:NUMtrials    % loop through trials
    inx = time > endtimes(iT);
    spk(iT,inx) = NaN;   % NaN trials after following event
end