function quickanalysis_gamblingtask(animalNO,sessionID,sessionspec)

%   Analysis of tetrode data, is designed as an offline analysis tool for
%   tetrode data and behavior, which can be executed in an unsupervised
%   manner on a daily bases. It gives a quick overview of the experiment
%   including response profiles of clustered neurons, light-triggered PSTH
%   and psychometric plots of behavior. It relies on CellBase data handling
%   system.
%
%   quickanalysis_gamblingtask(ANIMALNO,SESSIONID,SESSIONSPEC) performs the
%   analysis for a session specified by the first two input arguments.
%   SESSIONSPEC should be a 1x3 logical array indicating the need for
%   behavior, recording and stimulation analysis. Analyses the licking
%   behaviour of the animal and the neurons response to sound/ reward.
%
%   See also ADDNEWCELLS, PREALIGNSPIKES and VIEWCELL2B.
% -------------------------------------------------------------------------
% Cecília Pardo-Bellver, 2021
% Laboratory of Network Neurophysiology
% Institute of Experimental Medicine, Hungary.
% -------------------------------------------------------------------------

bigfig    = 0; % Make figure all screen
figformat = '.jpg'; % Change image file format
breath    = 0.01; % Pause time to visualise images
reTE      = 0; % Remake the TE files
% Input arguments ---------------------------------------------------------
narginchk(0,4)

% Animal, session, Behavior, recording or both
if nargin < 1
    disp('ERROR: No AnimalName provided.');
else
    if ischar(animalNO)
        animalID = animalNO;
    else
        animalID = ['BAd' num2str(animalNO,'%02.0f')];
    end
end
if nargin < 2
    disp('ERROR: No sessionID provided.');
end
if nargin < 3 
    isbeh  = 1;
    isrec  = 0;
    isstim = 0;
else
    isbeh  = sessionspec(1);
    isrec  = sessionspec(2);
    isstim = sessionspec(3);
end

fullpth = fullfile(getpref('cellbase','datapath'),animalID,sessionID);
cd(fullpth)

% Stop if error
dbstop if error

% Directories -------------------------------------------------------------
global DATAPATH %#ok<NUSED>
resdir = [getpref('cellbase','datapath') '\_response_profiles\' animalID '\'];
if ~isfolder(resdir)
    mkdir(resdir)
end
resdir2 = [getpref('cellbase','datapath') '\_behavior\' animalID '\'];
if ~isfolder(resdir2)
    mkdir(resdir2)
end

% Create Trial Events structure -------------------------------------------
if isbeh 
    if exist('TE_behaviour.mat','file') ~= 0 && reTE == 0
        TE_behaviour      = load('TE_behaviour.mat');
    elseif exist('TE_behaviour.mat','file') == 0 || reTE == 1
        behfile = dir([fullpth filesep animalID '_GamblingTask_'...
            sessionID(1:end-1) '*.mat']);
        behfile = behfile.name;
        TE_behaviour = TE_gamblingtask(fullfile(fullpth,behfile),1);
    end
   
    if isrec
        if exist('TE_recording.mat','file') ~= 0 && reTE == 0
             
        elseif exist('TE_recording.mat','file') == 0 || reTE == 1
            MakeTrialEvents_gambling(fullpth);% synchronize
        end   
    end
end

if isstim
    if exist('TE_StimEvent.mat','file') ~= 0 && reTE == 0
        
    elseif exist('TE_StimEvent.mat','file') == 0 || reTE == 1
        MakeTrialEvents_gambling(fullpth,'stim')
    end
end

% Update CellBase ---------------------------------------------------------
if isrec
    addnewcells('dir', [animalID '\' sessionID '\'], 'prealign', true);
    cellids = findcell('mouse',animalID,'session',sessionID);
    disp(cellids)
end

% Behavior ----------------------------------------------------------------
if isbeh
    StageNum = bynumbers_gamblingtask(animalID,sessionID,'s');
    stg = num2str(TE_behaviour.TrainingStage(1,1));
    switch (TE_behaviour.TrainingStage(1,1))	
     case 1
        % Lick raster 
	    H = figure;
        TE_behaviour.StimulusOn = zeros(1,length(TE_behaviour.TrainingStage));
        viewlickRL({animalID sessionID},'TE',TE_behaviour,'TriggerName',...
            'StimulusOn','SortEvent','TrialStart','eventtype','behav',...
            'ShowEvents',{{'StimulusOn'}},'Partitions','#TrialType',...
            'window',[-5 5],'LickSide','all');
	    fnm = fullfile(resdir2,[animalID '_' sessionID '_Stg' stg '_Licks' figformat]); 
        saveas(H,fnm)
        if bigfig == 1
            maximize_figure(H)
        end
        close(H)
        
        disp(['Current stage: Stg ' stg]);
        disp(['Advance new stage?[Y = 1, N = 0]: ' num2str(StageNum.NextStg)]);
    
     case 2
        % Lick raster all licks, pertition trial
	    H = figure;
	    viewlickRL({animalID sessionID},'TE',TE_behaviour,'TriggerName',...
            'StimulusOn','SortEvent','TrialStart','eventtype','behav',...
            'ShowEvents',{{'StimulusOn'}},'Partitions','#TrialType',...
            'window',[-5 5],'LickSide','all');
        fnm = fullfile(resdir2,[animalID '_' sessionID '_Stg' stg '_Licks' figformat]);
        saveas(H,fnm)
        if bigfig == 1
            maximize_figure(H)
        end
        close(H)
      
        % Lick raster R
	    H = figure;
	    viewlickRL({animalID sessionID},'TE',TE_behaviour,'TriggerName',...
            'StimulusOn','SortEvent','TrialStart','eventtype','behav',...
            'ShowEvents',{{'StimulusOn'}},'Partitions','#TrialType',...
            'window',[-5 5],'LickSide','R');
	    fnm = fullfile(resdir2,[animalID '_' sessionID '_Stg' stg '_LicksR' figformat]);
        saveas(H,fnm)
        if bigfig == 1
            maximize_figure(H)
        end
        close(H)
    
        % Lick raster L
	    H = figure;
	    viewlickRL({animalID sessionID},'TE',TE_behaviour,'TriggerName',...
            'StimulusOn','SortEvent','TrialStart','eventtype','behav',...
            'ShowEvents',{{'StimulusOn'}},'Partitions','#TrialType',...
            'window',[-5 5],'LickSide','L');
        fnm = fullfile(resdir2,[animalID '_' sessionID '_Stg' stg '_LicksL' figformat]);
        saveas(H,fnm)
        if bigfig == 1
            maximize_figure(H)
        end
        close(H)
        
        % Lick raster TrialT1
        H = figure;
	    viewlickTT({animalID sessionID},'TE',TE_behaviour,'TriggerName',...
            'StimulusOn','SortEvent','TrialStart','eventtype','behav',...
            'ShowEvents',{{'StimulusOn'}},'Partitions','#FirstLick',...
            'window',[-5 5],'LickSide','all','TrialType','1');
        fnm = fullfile(resdir2,[animalID '_' sessionID '_Stg' stg '_TrialT1' figformat]); 
        saveas(H,fnm)
        if bigfig == 1
            maximize_figure(H)
        end
        close(H)
        
        % Lick raster TrialT2
        H = figure;
	    viewlickTT({animalID sessionID},'TE',TE_behaviour,'TriggerName',...
            'StimulusOn','SortEvent','TrialStart','eventtype','behav',...
            'ShowEvents',{{'StimulusOn'}},'Partitions','#FirstLick',...
            'window',[-5 5],'LickSide','all','TrialType','2');
        fnm = fullfile(resdir2,[animalID '_' sessionID '_Stg' stg '_TrialT2' figformat]); 
        saveas(H,fnm)
        if bigfig == 1
            maximize_figure(H)
        end
        close(H)
        
        disp(['Current stage: Stg ' stg]);
        disp(['Advance new stage?[Y = 1, N = 0]: ' num2str(StageNum.NextStg)]);
        
     case 3       
        % Lick raster TrialT1
        H = figure;
	    viewlickTT({animalID sessionID},'TE',TE_behaviour,'TriggerName',...
            'StimulusOn','SortEvent','TrialStart','eventtype','behav',...
            'ShowEvents',{{'StimulusOn'}},'Partitions','#FirstLick',...
            'window',[-5 5],'LickSide','all','TrialType','1');
        fnm = fullfile(resdir2,[animalID '_' sessionID '_Stg' stg '_TrialT1' figformat]);
        saveas(H,fnm)
        if bigfig == 1
            maximize_figure(H)
        end
        close(H)
        
        % Lick raster TrialT2
        H = figure;
	    viewlickTT({animalID sessionID},'TE',TE_behaviour,'TriggerName',...
            'StimulusOn','SortEvent','TrialStart','eventtype','behav',...
            'ShowEvents',{{'StimulusOn'}},'Partitions','#FirstLick',...
            'window',[-5 5],'LickSide','all','TrialType','2');
	    fnm = fullfile(resdir2,[animalID '_' sessionID '_Stg' stg '_TrialT2' figformat]); 
        saveas(H,fnm)
        if bigfig == 1
            maximize_figure(H)
        end
        close(H)
        
    end   
end

% Neurons response profiles -----------------------------------------------
if isbeh && isrec   
    % Prealign spikes for trial events
    problem_behav_cellid = [];

    for iC = 1:length(cellids)
        cellid = cellids(iC);
        if exist(['EVENTSPIKES' cellid{1,1}(end-2) '_' cellid{1,1}(end) '.mat'],'file') == 0 ...
                || reTE == 1
            MakeTrialEvents_gambling(fullpth,'stim')
            try
                prealignSpikes_gambling(cellid,'FUNdefineEventsEpochs',...
                    @defineEventsEpochs_gambling,'filetype','event',...
                    'ifsave',1,'ifappend',0)
            catch
                disp('Error in prealignSpikes.');
                problem_behav_cellid = [problem_behav_cellid cellid];
            end
        end
        
    end
    
    % Aligning to StimulusOnset, Partition Reward 
    for k = 1:length(cellids)
        H = figure;
        viewcell2b_gambling(cellids(k),'TriggerName','StimulusOn','SortEvent',...
            'TrialStart','sigma', 0.07,'eventtype','behav','ShowEvents',...
            {{'DeliverAllFeedback'}},'Partitions','#Reward','window',[-3 3],...
            'PlotZeroLine','on')
        
        cellidt = cellids{k};
        cellidt(cellidt=='.') = '_';
        fnm = fullfile(resdir,[cellidt '_StimulusOn' figformat]);  
        saveas(H,fnm)
        pause(breath)
        if bigfig == 1
            maximize_figure(H)
        end
        close(H)
    end
    
    % Aligning to Feedback, Partition Reward
    for k = 1:length(cellids)
        H = figure;
        viewcell2b_gambling(cellids(k),'TriggerName','DeliverFeedback',...
            'SortEvent','TrialStart','sigma', 0.07,'eventtype','behav',...
            'ShowEvents',{{'StimulusOn'}},'Partitions','#Reward','window',[-3 3],...
            'PlotZeroLine','on')
        
       
        cellidt = cellids{k};
        cellidt(cellidt=='.') = '_';
        fnm = fullfile(resdir,[cellidt '_Reward' figformat]);   % save
        saveas(H,fnm)
        pause(breath)
        if bigfig == 1
            maximize_figure(H)
        end
        close(H)
    end
    
end

% Light effects -----------------------------------------------------------
if isrec && isstim
    % Prealign spikes to stimulus events
    problem_stim_cellid = [];
    for iC = 1:length(cellids)
        cellid = cellids(iC);
        if exist(['STIMSPIKES' cellid{1,1}(end-2) '_' cellid{1,1}(end) '.mat'],'file') == 0 ...
            || reTE == 1
            MakeTrialEvents_gambling(fullpth,'stim')
            try
                prealignSpikes_gambling(cellid,'FUNdefineEventsEpochs',...
                    @defineEventsEpochs_gambling,'filetype','stim',...
                    'ifsave',1,'ifappend',0)
            catch
                disp('Error in prealignSpikes.');
                problem_stim_cellid = [problem_stim_cellid cellid]; %#ok<AGROW>
            end
        end
    end
    
%     % View light-triggered raster and PSTH 
%     EventColors = [[0.4940 0.1840 0.5560];[0.9290 0.6940 0.1250];...
%         [0.6350 0.0780 0.1840]]; % Colors for PulseOn, PulseOff, BurstOff
    EventColors = [[1 0 0];[0 1 0];[0 0 1]]; 
    EventColors = mat2cell(EventColors,ones(size(EventColors,1),1),3);
    for iC = 1:length(cellids)
        cellid = cellids(iC);
        H = figure;
        viewcell2b_gambling(cellid,'TriggerName','BurstOn','SortEvent','BurstOff',...
            'ShowEvents',{{'PulseOn','PulseOff','BurstOff'}},...
            'ShowEventsColors',{EventColors},'FigureNum',H,...
            'eventtype','stim','window',[-0.2 0.5],'dt',0.001,'sigma',0.001,...
            'PSTHstd','on','Partitions','#PulseID',...
            'EventMarkerWidth',0,'PlotZeroLine','off');
        
        cellidt = cellid{1};
        cellidt(cellidt=='.') = '_';
        fnm = fullfile(resdir,[cellidt '_PSTH' figformat]); 
        saveas(H,fnm)
        pause(breath)
        if bigfig == 1
            maximize_figure(H)
        end
        close(H)
    end
end

% % Cluster quality
% % global CheckCluster_DrawAllWaves
% if isrec
%     CheckCluster_DrawAllWaves = true; %#ok<NASGU>
%     BatchSessionClust(fullpth)
%     CheckCluster_DrawAllWaves = false;
% end

end