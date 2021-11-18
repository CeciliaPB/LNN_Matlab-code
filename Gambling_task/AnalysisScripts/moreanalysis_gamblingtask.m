function moreanalysis_gamblingtask(animalNO,sessionID,sessionspec)

%   Analysis of tetrode data, is designed as an offline analysis tool for
%   tetrode data and behavior, which can be executed in an unsupervised
%   manner on a daily bases. It gives a further overview of the experiment
%   including response profiles of clustered neurons, light-triggered PSTH
%   and psychometric plots of behavior. It relies on CellBase data handling
%   system. Complements the script analysis_gamblingtask, that analyses the
%   licking behaviour of the animal and the neurons response to sound/
%   reward.
%
%   analysis_gamblingtask(ANIMALNO,SESSIONID,SESSIONSPEC) performs the
%   analysis for a session specified by the first two input arguments.
%   SESSIONSPEC should be a 1x3 logical array indicating the need for
%   behavior, recording and stimulation analysis.
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
    isrec  = 1;
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
            TE_recording = load('TE_recording.mat');
        elseif exist('TE_recording.mat','file') == 0 || reTE == 1
            MakeTrialEvents_gambling(fullpth);% synchronize
        end   
    end
end

% Update CellBase ---------------------------------------------------------
if isrec
    addnewcells('dir', [animalID '\' sessionID '\'], 'prealign', true);
    cellids = findcell('mouse',animalID,'session',sessionID);
    disp(cellids)
end

% Neurons response profiles -----------------------------------------------
if isbeh && isrec   
    % Prealign spikes for trial events
    problem_behav_cellid = [];
    for iC = 1:length(cellids)
        cellid = cellids(iC);
        if exist(['EVENTSPIKES' cellid{1,1}(end-2) '_' cellid{1,1}(end) '.mat'],'file') == 0 || reTE == 1
            MakeTrialEvents_gambling(fullpth,'stim')
            try
                prealignSpikes_gambling(cellid,'FUNdefineEventsEpochs',@defineEventsEpochs_gambling,...
                    'filetype','event','ifsave',1,'ifappend',0)
            catch
                disp('Error in prealignSpikes.');
                problem_behav_cellid = [problem_behav_cellid cellid];
            end
        end
    end
    
    % Aligning to StimulusOffset, Partition Reward 
    for k = 1:length(cellids)
        H = figure;
        viewcell2b_gambling(cellids(k),'TriggerName','StimulusOff','SortEvent','TrialStart',...
            'sigma', 0.07,'eventtype','behav','ShowEvents',{{'DeliverAllFeedback'}},...
            'Partitions','#Reward','window',[-3 3],'PlotZeroLine','on')
        
        cellidt = cellids{k};
        cellidt(cellidt=='.') = '_';
        fnm = fullfile(resdir,[cellidt '_StimulusOff' figformat]);   % save
        saveas(H,fnm)
        pause(breath)
        if bigfig == 1
            maximize_figure(H)
        end
        close(H)
    end
    
    % Aligning to StimulusOnset, Partition Punishment
    stg = num2str(TE_behaviour.TrainingStage(1,1));
    if stg == 3
        for k = 1:length(cellids)
            H = figure;
            viewcell2b_gambling(cellids(k),'TriggerName','StimulusOn',...
                'SortEvent','TrialStart','sigma', 0.07,'eventtype','behav',...
                'ShowEvents',{{'StimulusOn'}},'Partitions','#Punishment',...
                'window',[-3 3],'PlotZeroLine','on');
            
            cellidt = cellids{k};
            cellidt(cellidt=='.') = '_';
            fnm = fullfile(resdir,[cellidt '_Punishment' figformat]);   % save
            saveas(H,fnm)
            pause(breath)
            if bigfig == 1
                maximize_figure(H)
            end
            close(H)
        end
    end
    
    % Aligning to Feedback, Partition TrialType
    for k = 1:length(cellids)
        H = figure;
        viewcell2b_gambling(cellids(k),'TriggerName','DeliverAllFeedback',...
            'SortEvent','TrialStart','sigma', 0.07,'eventtype','behav',...
            'ShowEvents',{{'DeliverAllFeedback'}},'Partitions','#TrialType',...
            'window',[-3 3],'PlotZeroLine','on');
        
        cellidt = cellids{k};
        cellidt(cellidt=='.') = '_';
        fnm = fullfile(resdir,[cellidt '_TrialType' figformat]);   % save
        saveas(H,fnm)
        pause(breath)
        if bigfig == 1
            maximize_figure(H)
        end
        close(H)
    end
    
end

% % Cluster quality
% global CheckCluster_DrawAllWaves
% if isrec
%     CheckCluster_DrawAllWaves = true;
%     BatchSessionClust(fullpth)
%     CheckCluster_DrawAllWaves = false;
% end

end