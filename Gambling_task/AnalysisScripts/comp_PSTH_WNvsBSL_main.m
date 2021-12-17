function comp_PSTH_WNvsBSL_main(resdir,varargin)

% Checks if a cell changes significantly firing rate between a certain time
% window (the test window) and a baseline period. 
% NOTE: we do not assume that baseline and test windows are adjacent and if
% in the test window both excitation and inhibition happens, we consider the
% longer of the two.
%
%   Input parameters:
%    - resdir, where to save plots and graphs
%    - varargin, See Default arguments.
% It can be used with or without CellBase, default NO.
%
% See also comp_PSTH_WNvsBSL_statszero, comp_PSTH_WNvsBSL_stats and
% comp_PSTH_WNvsBSL_stats_CB.
% -------------------------------------------------------------------------
% Based on PSTH_stat_test_vs_bsln_length by Nicola Solari
% 
% Cecília Pardo-Bellver, 2021
% Laboratory of Network Neurophysiology
% Institute of Experimental Medicine, Hungary.
% -------------------------------------------------------------------------

% Default arguments
prs = inputParser;
addRequired(prs,'resdir',@(s)isempty(s)|isfolder(s)) % Directory

addParameter(prs,'CellBase','false') % Default NOT use CellBase. 
addParameter(prs,'CellList',[]) % CellList to load, a .xls file.
addParameter(prs,'TTL',[]) % Variable containing the TTLs, a .mat file.
addParameter(prs,'SessionID','') % Default all sessions. 
addParameter(prs,'AnimalID','') % Default all animals. 
addParameter(prs,'issave',0) % Default not save. To save change to 1
addParameter(prs,'responsedir',{}) % Results folder
addParameter(prs,'Quality',0) % Use Quality to select cells
addParameter(prs,'partitions',{'#Reward', '#Punishment'}) % Partitions:'#TrialType' Whatever you want/have
addParameter(prs,'PoI',{'Sound','Shock','Laser'}) % Periods of interest: these are names, the specs are after
addParameter(prs,'sig_thr',0.001,@isnumeric) % Significance threshold you want to save in TheMatrix
parse(prs,resdir,varargin{:})
g = prs.Results;

issave  = g.issave;
sig_thr = g.sig_thr;

switch g.CellBase
    case 'false'
        for poi = 1 : length(g.PoI) % PoI loop
            % Choose CellBase & Cells
            if isempty(g.CellList)
                error('No CellList provided!')
            end
            Tab = readtable(g.CellList);
            CellList = Tab;
            if isempty(g.responsedir)
                responsedir = fullfile(resdir,'WNvsBSL'); % Where to save stats and figs
            else
                responsedir = fullfile(resdir,g.responsedir);
            end
            if g.Quality ~= 0
                qualityids = (Tab.Quality< g.Quality); % Select cells based on quality
                Tab = Tab(qualityids,:);
            else
                qualityids = (Tab.Lratio<0.15 & Tab.IsolationDistance>20); % Select cells with good lratio and ID (Sergio knows)
                Tab = Tab(qualityids,:);
            end
            if ~isempty(g.AnimalID)
                animalids = (contains(Tab.AnimalID, g.AnimalID) == 1);
                Tab = Tab(animalids,:);
            end
            if ~isempty(g.SessionID)
                sessionids = Tab.CellID(contains(Tab.SessionID, g.SessionID) == 1);
                Tab = Tab(sessionids,:);
            end
            if ~isfolder(responsedir)
                mkdir(responsedir)
            end
            cellids = Tab.CellID;
            
            % PSTH general parameters
            wn    = [-3 3]; % Raster window in seconds
            dt    = 0.001; % Raster resolution, in seconds
            sigma = 0.02; % Controls smoothing for 'spsth'
            focus = g.PoI{poi};
            
            switch focus % Define the "periods of interest" according to my protocol
                case 'Sound'
                    % Comparison to sound delivery
                    BLevent = 'Sound'; % Event used for baseline BL
                    bwin  = [-2 0]; % BL window for MW-test: the 2 sec before the BLevent
                    twin  = [0 1]; % WN for MW-test:  the 1 sec after the BLevent
                    if sum(abs(bwin+twin))> 1
                        sigma = sigma/4;
                    end
                case 'Shock'
                    % Comparison to shock delivery
                    BLevent = 'Shock';
                    bwin  = [-0.4 0];
                    twin  = [0 0.2];
                case 'Laser'
                    % Comparison to laser delivery
                    BLevent = 'Laser';
                    bwin  = [-0.4 0];
                    twin  = [0 0.2];
            end
            
            % Add property for grouping
            propname = ['statsPSTH_' focus];
            CellList.(['CompID_' BLevent]) = cell(length(CellList.CellID),1);
            CellList.(['iMW_pv_' BLevent]) = nan(length(CellList.CellID),1);
            CellList.(['aMW_pv_' BLevent]) = nan(length(CellList.CellID),1);
            
            for iC = 1:length(cellids) % Cell loop
                cellid = cellids{iC}; % Current cell
                A = [Tab.AnimalID{iC},Tab.SessionID{iC}];
                B = length(A);
                cellTT = cellid(B+2:end);

                mainFolder = 'L:\_Neurons';
                folder     = [Tab.AnimalID{iC,:},Tab.SessionID{iC,:}];
                cd([mainFolder filesep folder]);
                
                clear events
                load(g.TTL,'events');
                event  = events.(BLevent);
                
                try
                    stats = rasterPSTH(cellTT,event,bwin,twin,wn,dt,sigma,...
                        issave,responsedir,propname,sig_thr,BLevent,cellid);
                    
                    ncell = find(ismember(CellList.CellID,cellid));
                    CellList.(['CompID_' BLevent])(ncell,:) = {stats.tag};
                    CellList.(['iMW_pv_' BLevent])(ncell,:) = stats.MWp_i;
                    CellList.(['aMW_pv_' BLevent])(ncell,:) = stats.MWp_a;
                                        
                catch
                    disp (['Error working on ' cellid])
                    close all
                end
            end
            cd(mainFolder);
            if issave
                writetable(CellList,g.CellList)
            end
        end % PoI loop
        
    case 'true'
        for part = 1 : length(g.partitions) % Partitions loop
            for poi = 1 : length(g.PoI) % PoI loop
                % Choose CellBase & Cells
                choosecb('Risk_ass_base');
                loadcb
                responsedir = fullfile(resdir,'WNvsBSL'); % Where to save stats and figs
                if g.Quality ~= 0
                    cellids = selectcell(['"Quality"<' g.Quality]); % Select cells based on quality
                else
                    cellids = selectcell('"ID">20&"Lratio"<0.15'); % Select cells with good lratio and ID (Sergio knows)
                end
                if ~isempty(g.AnimalID)
                    cellids = cellids(contains(cellids, g.AnimalID) == 1);
                end
                if ~isempty(g.SessionID)
                    cellids = cellids(contains(cellids, g.SessionID) == 1);
                end
                if ~isfolder(responsedir)
                    mkdir(responsedir)
                end
                
                % PSTH general parameters
                wn    = [-3 3]; % Raster window in seconds
                dt    = 0.001; % Raster resolution, in seconds
                sigma = 0.02; % Controls smoothing for 'spsth'
                focus = g.PoI{poi};
                
                switch focus % Define the "periods of interest" according to my protocol
                    case 'Sound'
                        % Comparison to cue delivery
                        BLevent = 'StimulusOn'; % Event used for baseline BL
                        bwin    = [-0.4 0]; % BL window for MW-test: the 0.4 sec before the BLevent
                        WNevent = 'StimulusOn'; % Event used for window WN
                        twin    = [0 0.2]; % WN for MW-test:  the 0.2 sec after the WNevent
                        event   = {BLevent,WNevent};
                        partition = g.partitions{part};   % partition trials
                    case 'Feedback'
                        % Comparison to reward delivery
                        BLevent = 'StimulusOn';
                        bwin    = [-0.4 0];
                        WNevent = 'DeliverAllFeedback';
                        twin    = [0 0.2];
                        event   = {BLevent,WNevent};
                        partition = g.partitions{part};
                    case 'Feedback2'
                        % Comparison to reward delivery, BL just before
                        BLevent = 'DeliverAllFeedback';
                        bwin    = [-0.2 0];
                        WNevent = 'DeliverAllFeedback';
                        twin    = [0 0.2];
                        event   = {BLevent,WNevent};
                        partition = g.partitions{part};
                end
                
                % Add property for grouping
                propname = ['statsPSTH_' focus '_' partition(2:end)]; % property name saved in ANALYSES
                if ~ismember(propname,listtag('prop'))
                    insertdata([CELLIDLIST' num2cell(nan(size(CELLIDLIST')))],...
                        'type','property','name',propname)
                end
                clear xlsdata
                xlsdata = table;
                
                for iC = 1:length(cellids) % Cell loop
                    cellid = cellids{iC}; % Current cell
                    try
                        xlsdata.CellID(iC,:)   = cellid;
                        xlsdata.FociPart(iC,:) = propname;
                        
                        stats = rasterPSTH_CB(cellid,event,bwin,twin,partition,...
                            wn,dt,sigma,issave,responsedir,propname,sig_thr);
                        
                        %     Add activity change to CellBase: 1 = inhibited,
                        %     2 = neutral, 3 = excited, 9 = something is wrong
                        a = nan(1,length(stats));
                        for nst = 1:length(stats) % Check all partitions#
                            if     stats{nst}.MWp_i > sig_thr && stats{nst}.MWp_a < sig_thr
                                a(nst) = 3; % activated
                            elseif stats{nst}.MWp_i < sig_thr && stats{nst}.MWp_a > sig_thr
                                a(nst) = 1; % inhibited
                            elseif stats{nst}.MWp_i > sig_thr && stats{nst}.MWp_a > sig_thr
                                a(nst) = 2; % nada
                            elseif stats{nst}.MWp_i < sig_thr && stats{nst}.MWp_a < sig_thr
                                if abs(stats{nst}.inhibition_end) - abs(stats{nst}.inhibition_start)...
                                        < abs(stats{nst}.activation_end) - abs(stats{nst}.activation_start) %if longer activated
                                    a(nst) = 3;   % activated
                                else
                                    a(nst) = 1;   % inhibited
                                end
                            else
                                a(nst) = 9;
                            end
                            a(nst) = a(nst)*10^(length(stats)-nst);
                        end
                        
                        if ~isempty(find((a == 0),1)) % something went wrong
                            disp (['Error using setvalue for ' cellid])
                        else
                            actchange = sum(a);
                            setvalue(cellid,propname,actchange);
                        end
                        
                        for nst = 1:length(stats) % Check all partitions#
                            xlsdata.(['CompID' num2str(nst)])(iC,:) = stats{nst}.tag;
                            xlsdata.(['iMW_pv' num2str(nst)])(iC,:) = stats{nst}.MWp_i;
                            xlsdata.(['aMW_pv' num2str(nst)])(iC,:) = stats{nst}.MWp_a;
                        end
                        
                    catch
                        disp (['Error working on ' cellid])
                        close all
                    end
                end
                if issave
                    fname = (['WNvsBSL_data_',propname,'.xls']);
                    writetable(xlsdata,fname)
                end
            end % PoI loop
        end % Partition loop
end

% -------------------------------------------------------------------------
function stats = rasterPSTH(cellTT,event,bwin,twin,wn,dt,sigma,issave,responsedir,propname,sig_thr,tag,cellid)

% Raster plot baseline and window of interest
figure(1) % Baseline
raster_plot(cellTT,event,'eventline',1,'window',bwin)

% psth = nanmean(binraster) / dt;
% fs       = 30000; % sampling frequency
% psth_bin = 600; % 600 = 20ms
% psth = ttl_psth (stimes*fs, event*fs, psth_bin, 'pre', abs(wn(1)),'post', abs(wn(2)));
% [psth_spx, psth_t] = psth_hist(psth, psth_bin);

figure(2)
raster_plot(cellTT,event,'eventline',1,'window',twin)

% Peri-event time histogram and stats generation
stats = comp_PSTH_WNvsBSL_stats(cellTT,event,bwin,twin,'dt',dt,'sigma',sigma,...
    'relative_threshold',0.1,'tag',tag,'sig_thr',sig_thr,'window',wn); % calculate psth

figure(5)
f1   = figure(1);
h1   = get(f1,'Children');
wdif = (sum(abs(bwin))+sum(abs(twin)))/sum(abs(twin));
wdif2 = sum(abs(bwin))/sum(abs(twin));
newh = copyobj(h1,5); % partition nn scatter
ax   = get(figure(5),'Children');
ax(1).XColor = 'none';
pos1 = get(newh,'Position');
set(newh,'Position',[pos1(1),pos1(2)*5,pos1(3)*wdif2/wdif,pos1(4)/2]);
title([cellid ': ' propname],'Interpreter','none','HorizontalAlignment','left');

f2   = figure(2);
h2   = get(f2,'Children');
newh = copyobj(h2,5); % partition nn scatter
ax   = get(figure(5),'Children');
ax(1).YColor = 'none';
ax(1).XColor = 'none';
pos2 = get(newh,'Position');
set(newh,'Position',[pos2(1)+(pos1(3)*wdif2/wdif),pos2(2)*5,pos2(3)/wdif,pos2(4)/2]);

f3   = figure(3);
h3    = get(f3,'Children');
newh = copyobj([h3(1) h3(2)],5);
figure(5)
set(newh(2),'Position',[pos2(1),pos2(2),pos2(3),pos2(4)/2]);

figure(5)
maximize_figure
pause(0.01)

% Save figure and stats
if issave
    fnm     = fullfile(responsedir,[cellid '_' propname '.jpg']);
    saveas(figure(5),fnm)
    fnm2    = fullfile(responsedir,[cellid '_' propname '.mat']);
    save(fnm2,'stats')
end

close all

% -------------------------------------------------------------------------
function stats = rasterPSTH_CB(cellid,event,bwin,twin,partition,wn,dt,sigma,issave,responsedir,propname,sig_thr)

cellpath = 'L:\Cecilia\Behaviour\Risk assesment';
% Raster plot and PSTH
try
    TE = load([cellpath '\' cellid(1:5) '\' cellid(7:15) '\' 'TE_recording.mat']); % load trial events
catch
    disp (['No recording on ' cellid])
end
SP = loadcb(cellid,'EVENTSPIKES');
fld = fieldnames(TE);
if ~isequal(length(SP.event_stimes{1}),length(TE.(fld{1})))
    error('MATLAB:trialMismatch',...
        'Trial number mismatch between TrialEvents and EVENTSPIKES.')
end

% Raster plot baseline and window of interest
figure(1)
viewcell2b_gambling(cellid,'TriggerName',event{1},'SortEvent','TrialStart',...
    'sigma',sigma,'eventtype','behav','ShowEvents',event{1},...
    'Partitions',partition,'window',bwin, 'PSTHPlot',false,'PrintCellID','off'); % 'PSTHstd', 'on',
figure(2)
viewcell2b_gambling(cellid,'TriggerName',event{2},'SortEvent','TrialStart',...
    'sigma',sigma,'eventtype','behav','ShowEvents',event{2},...
    'Partitions',partition,'window',twin, 'PSTHPlot',false,'PrintCellID','off'); % 'PSTHstd', 'on',

% Peri-event time histogram and stats generation
stats = comp_PSTH_WNvsBSL_stats_CB(cellid,event,bwin,twin,'dt',dt,'sigma',sigma,...
    'partitions',partition,'isadaptive',0,'maxtrialno',Inf,'relative_threshold',...
    0.1,'display',true,'sig_thr',sig_thr,'window',wn); % calculate psth

figure(5)
cellidt = regexprep(cellid,'\.','_');
f1   = figure(1);
h1   = get(f1,'Children');
nfig = length(h1);
wdif = (sum(abs(bwin))+sum(abs(twin)))/sum(abs(twin));
wdif2 = sum(abs(bwin))/sum(abs(twin));
for nn = 1:floor(nfig/2)
    copyobj(h1(nn),5); % partition nn color 
    newh = copyobj(h1(floor(nfig/2)+nn+1),5); % partition nn scatter
    ax   = get(figure(5),'Children');
    ax(1).YColor = 'none';
    pos1 = get(newh,'Position');
    set(newh,'Position',[pos1(1),pos1(2),pos1(3)*wdif2/wdif,pos1(4)]);
    if nn == floor(nfig/2) 
        title([cellidt ': ' propname],'Interpreter','none','HorizontalAlignment','left');
    end
end

f2   = figure(2);
h2   = get(f2,'Children');
nfig = length(h2);
for nn = 1:floor(nfig/2)
    copyobj(h2(nn),5); % partition nn color 
    newh = copyobj(h2(floor(nfig/2)+nn+1),5); % partition nn scatter
    pos2 = get(newh,'Position');
    set(newh,'Position',[pos2(1)+(pos1(3)*wdif2/wdif),pos2(2),pos2(3)/wdif,pos2(4)]);
end

f3   = figure(3);
h3    = get(f3,'Children');
newh = copyobj([h3(1) h3(2)],5); % part 1 scatter
figure(5)
set(newh(2),'Position',[pos2(1),pos2(2)/(nfig*2),pos2(3),pos2(4)]);

figure(5)
maximize_figure
pause(0.01)

% Save figure and stats
if issave
    fnm     = fullfile(responsedir,[cellidt '_' propname '.jpg']);
    saveas(figure(5),fnm)
    fnm2    = fullfile(responsedir,[cellidt '_' propname '.mat']);
    save(fnm2,'stats')
end

close all