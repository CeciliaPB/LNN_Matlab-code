function viewlickTT(cbID,varargin)
%VIEWLICK   Lick raster and PSTH.
%   VIEWLICK(CBID,'TRIGGERNAME',TRIGEVENT,'SORTEVENT',SEVENT,'EVENTTYPE
%   ',EVTYPE) plots an event triggered lick raster and PSTH for a given
%   session (CBID: 1-by-2 cell containing animal ID and session ID).
%   TRIGEVENT is used as trigger and trials are sorted according to SEVENT.
%   'EVENTTYPE' can be 'STIM' or 'BEHAV'. Optional input argument pairs:
%       'ShowEvets', events indicated on the raster plots
%       'ShowEventsColors', specifies colors for the events shown
%       'FigureNum', specifies the figure handle to plot on
%       'window', window size for the plots
%       'dt', time resoultion
%       'sigma', determines the smoothing kernel for the smoothed PSTH (see
%           SMOOTHED_PSTH)
%       'isadaptive', default = false - if true, adaptive PSTH algorithm is
%           used (see APSTH and BINRASTER2APSTH)
%       'PSTHstd', 'on' or 'off', controls whether SD is plotted on the PSTH
%           panel
%       'Partitions', e.g. 'all', '#ResponseType', can be used to partition
%           the trials according to different variables, multiple rasters 
%           are plotted
%       'EventMarkerWidth', specifies marker size for events shown
%           'PlotZeroLine', 'on' or 'off', controls whether a line
%           indicating zero appears on the plots
%
%   For example calls, see QUICKANALISYS2.
%
%   See also VIEWCELL2B.
% -------------------------------------------------------------------------
%   Cec�lia Pardo-Bellver, 2021
%   Laboratory of Network Neurophysiology
%   Institute of Experimental Medicine, Hungary.
%
%   Based on VIEWCELL by Balazs Hangya   Edit log: BH 5/2/14
% -------------------------------------------------------------------------

% Default arguments
default_args={...
    'window',               [-0.5 1];...
    'dt',                   0.01;...
    'sigma',                0.02;...
    'isadaptive'            false;...
    'FigureNum',            1;...
    'TriggerName',          'PulseOn';...
    'SortEvent',            'PulseOn';...
    'eventtype',            'stim';... % 'behav'
    'ShowEvents',           {{'PulseOn'}};...
    'ShowEventsColors',     {{[0 0.8 0] [0.8 0.8 0] [0 0.8 0.8] [0.8 0 0.8] [0.8 0 0] [0 0 0.8]}};...
    'Num2Plot',             'all';...
    'PlotDashedEvent',      '';...
    'PlotDashedCondition',  'min';...
    'PSTHPlot',             1;...
    'PSTHlinewidth',        1.5;...
    'DashedLineStyle',      ':';...
    'LastEvents',           '';...
    'Partitions',           'all';...
    'PrintCellID',          'on';...
    'PrintCellIDPos',       'bottom-right';...
    'BurstPSTH',            'off';...
    'LickSide',             'all';... % Choose 'all', 'R' right, 'L' left
    'TE',                   '';... % Default not loaded
    'TrialType',            'all';... % Chose 'all', '1' TrialType 1, '2' TT2
    };
[g,error] = parse_args(default_args,varargin{:}); %#ok<ASGLU>

% Load trial events (unsynchronized)
animalID = cbID{1,1};
sessionID = cbID{1,2};
cbdir = getpref('cellbase','datapath');
switch g.eventtype
    case 'stim'
        datapath = fullfile(cbdir,animalID,sessionID,'StimEvents.mat');
        if isempty(g.TE)
            TE = load(datapath);
        elseif ~isempty(g.TE)
            TE = g.TE;
        end
    case {'event','behav'}
        datapath = fullfile(cbdir,animalID,sessionID,'TE_behaviour.mat');
        if isempty(g.TE)
            TE = load(datapath);
        elseif ~isempty(g.TE)
            TE = g.TE;
        end
end

%--------------------------------------------------------------------------
% Preprocessing
%--------------------------------------------------------------------------

% Time
margin = g.sigma * 3;     % add an extra margin to the windows
time   = g.window(1)-margin:g.dt:g.window(2)+margin;   % time base array

% Lick times, depends on g.LickSide ---------------------------------------
switch (g.LickSide)
    case 'all' 
        NumTrials = length(TE.LickIn);
        LickIn = cell(1,NumTrials);
        for ii = 1:NumTrials
            if TE.FirstLick(1,ii) == 1
                LickIn(1,ii) = TE.LickRIn(1,ii);
            elseif TE.FirstLick(1,ii) == 2
                LickIn(1,ii) = TE.LickLIn(1,ii);
            end
        end
        stimes = arrayfun(@(s)LickIn{s}-TE.(g.TriggerName)(s),1:NumTrials,...
            'UniformOutput',false);
    case 'R'
        NumTrials = length(TE.LickRIn);
        stimes = arrayfun(@(s)TE.LickRIn{s}-TE.(g.TriggerName)(s),1:NumTrials,...
            'UniformOutput',false);
    case 'L'
        NumTrials = length(TE.LickLIn);
        stimes = arrayfun(@(s)TE.LickLIn{s}-TE.(g.TriggerName)(s),1:NumTrials,...
            'UniformOutput',false);
end

switch (g.TrialType)
    case 'all' 
        % No change in stimes
    case '1'
        TT1 = cell(1,NumTrials);
        for kk = 1:NumTrials
            if TE.TrialType(kk) == 1
                TT1(1,kk) = stimes(1,kk);
            end
        end
        stimes = TT1;
    case '2'
        TT2 = cell(1,NumTrials);
        for kk = 1:NumTrials
            if TE.TrialType(kk) == 2
                TT2(1,kk) = stimes(1,kk);
            end
        end
        stimes = TT2;
end

% Event windows
window_margin = g.window;
if iscell(TE.(g.TriggerName))
    ev_windows = cellfun(@(s)repmat(g.window',1,size(s,1)),TE.(g.TriggerName),...
        'UniformOutput',false);
else
    ev_windows = repmat(g.window',1,NumTrials);
end

%--------------------------------------------------------------------------
% Make the main raster
%--------------------------------------------------------------------------

% Calculate binraster
NUMtrials = length(TE.(g.SortEvent));
if iscell(stimes{1})   % deal with lick-aligned raster
    stimes2 = [stimes{1:end}];
    binraster0 = stimes2binraster(stimes2,time,g.dt);
    binraster = nan(NUMtrials,size(binraster0,2));
%     binraster2 = nan(NUMtrials,size(binraster0,2));
    for k = 1:NUMtrials   % calculate sum of rows for each trial, which will be used for the PSTH
        sind = sum(cellfun(@length,stimes(1:k-1))) + 1;
        eind = sind + length(stimes{k}) - 1;
%         disp([sind eind])
        binraster(k,:) = mean(binraster0(sind:eind,:),1);
%         binraster2(k,:) = sum(stimes2binraster(stimes{k},time,g.dt),1);
    end
else
    binraster = stimes2binraster(stimes,time,g.dt);
end

% For variable windows, change padding to NaN to ensure correct averaging - BH
if ~isempty(g.LastEvents)
    for iT = 1:NUMtrials    % loop through trials
        inx = time > ev_windows(iT,2);
        binraster(iT,inx) = NaN;
    end
end

% Partition trials
[COMPTRIALS, TAGS] = partition_trials(TE,g.Partitions);
vinx = cellfun(@(s)(~isempty(s)),COMPTRIALS);
COMPTRIALS = COMPTRIALS(vinx);
TAGS = TAGS(vinx);
trigev = TE.(g.TriggerName);
if ~iscell(trigev)
    valid_trials = find(~isnan(trigev));
else
    valid_trials = find(cellfun(@(s)~isempty(s),trigev));
end

% Calculate PSTH
switch g.isadaptive
    case {false,0}
        [psth, spsth, spsth_se] = binraster2psth(binraster,g.dt,g.sigma,...
            COMPTRIALS,valid_trials);
    case {true, 1}
        [psth, spsth, spsth_se] = binraster2apsth(binraster,g.dt,g.sigma,...
            COMPTRIALS,valid_trials);   % adaptive PSTH
    case 2
        [psth, spsth, spsth_se] = binraster2dapsth(binraster,g.dt,g.sigma,...
            COMPTRIALS,valid_trials);   % adaptive PSTH
end
EventTimes = trialevents2relativetime(TE,g.TriggerName,g.ShowEvents);

% Sort trials
NUMevents = length(g.SortEvent);
if iscellstr(g.SortEvent)
    sort_var = nan(NUMevents,NUMtrials);
    for iS = 1:NUMevents
        sort_var(iS,:) = TE.(g.SortEvent{iS}) - TE.(g.TriggerName);
    end
    sort_var = min(sort_var);
elseif ~isempty(g.SortEvent)
    if ~iscell(TE.(g.TriggerName))
        sort_var = TE.(g.SortEvent) - TE.(g.TriggerName);
    else
        gte = nan(1,NUMtrials);
        inx = ~cellfun(@isempty,TE.(g.TriggerName));
        gte(inx) = cell2mat(cellfun(@(s)s(1),TE.(g.TriggerName)(inx),...
            'UniformOutput',false));
        sort_var = TE.(g.SortEvent) - gte;
    end
else
    sort_var = NaN;
end

switch (g.LickSide)
    case 'all' 
        newTAGS = {'TrialType=1','TrialType=2'};
    case 'R'
        newTAGS = {'R_TrialType=1','R_TrialType=2'};
    case 'L'
        newTAGS = {'L_TrialType=1','L_TrialType=2'};
end

switch (g.TrialType)
    case 'all' 
        newTAGS = {'TrialType=1','TrialType=2'};
    case '1'
        newTAGS = {'TT1_R','TT1_L'};
    case '2'
        newTAGS = {'TT2_R','TT2_L'};
end

[mylabels, mycolors, mycolors2,mylinestyle] = makeColorsLabels(@defineLabelsColors_Cecilia,...
    newTAGS);
XLabel = ['Time - ' g.TriggerName];
YLabel = 'Rate (Hz)';

%--------------------------------------------------------------------------
% Raster + PSTH
%--------------------------------------------------------------------------

% Plot the raster
fhandle0 = plot_raster2a(stimes,time,valid_trials,COMPTRIALS,mylabels,...
    EventTimes,window_margin,ev_windows,sort_var,g,'Colors',{mycolors},...
    'Colors2',{mycolors2},'NumTrials2Plot',g.Num2Plot);
if isfield(g,'Legend')
    mylabels = g.Legend;
end

% Plot the PSTH
if g.PSTHPlot == 1
    if ~isempty(g.PlotDashedEvent)
        if nansum(TE.TrialStart) ~= 0   % i.e. TrialStart is an actual timestamp
            g.PlotDashedTime = eval(['TE.' g.PlotDashedEvent]);
        else   % if TrialStart is a dummy variable set to 0
            g.PlotDashedTime = ev_windows(2,valid_trials);
        end
        switch g.PlotDashedCondition
            case 'min'
                g.PlotDashedTime = nanmin(g.PlotDashedTime);
            case 'max'
                g.PlotDashedTime = nanmax(ev_windows(2,valid_trials));
            case 'mean'
                g.PlotDashedTime = nanmean(ev_windows(2,valid_trials));
            case 'median'
                g.PlotDashedTime = nanmedian(ev_windows(2,valid_trials));
        end
    end
    plot_timecourse(time,spsth,spsth_se,g,'FigureNum',fhandle0(end),...
        'Colors',{mycolors},'LineStyle',{mylinestyle},'Legend',{mylabels},...
        'XLabel',XLabel,'YLabel',YLabel);
    axis tight
end
if strcmpi(g.PrintCellID,'on')
    fstamp(sessionID,80,g.PrintCellIDPos);
end

% Link axes
A = findobj(allchild(gcf),'Type','axes');
Am = findobj(A,'YLim',[0 1]);
linkaxes(setdiff(A,Am),'x');