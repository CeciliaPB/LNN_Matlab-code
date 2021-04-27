function quickanalysis_gamblingtask(animalNO,sessionID,sessionspec,protocoltag)

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
%   behavior, recording and stimulation analysis.
%
%   quickanalysis_gamblingtask(ANIMALNO,SESSIONID,SESSIONSPEC,PROTOCOLTAG)
%   accepts a PROTOCOLTAG argument to allow calls to trial event conversion
%   programs for different versions of the behavioral protocol.
%
%   See also ADDNEWCELLS, PREALIGNSPIKES and VIEWCELL2B.
% -------------------------------------------------------------------------
% Based on quickanalysis_pavlovian2_p.m
% By Heguedus Panna.
%
% Modified by Cecília Pardo-Bellver, 2020
% Laboratory of Network Neurophysiology
% Institute of Experimantal Medicine, Hungary.
% -------------------------------------------------------------------------

% Input argument check
narginchk(0,4)

% Animal, session
if nargin < 1
    disp('ERROR: No AnimalName provided.');
else
    if ischar(animalNO)
        animalID = animalNO;
    else
        animalID = ['HB' num2str(animalNO)];
    end
end
if nargin < 2
    disp('ERROR: No sessionID provided.');
end

% Behavior, recording or both
if nargin < 4
    protocoltag = '';
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

% Directories
global DATAPATH %#ok<NUSED>
resdir = [getpref('cellbase','datapath') '\_response_profiles\' animalID '\'];
if ~isfolder(resdir)
    mkdir(resdir)
end
resdir2 = [getpref('cellbase','datapath') '\_behavior\' animalID '\'];
if ~isfolder(resdir2)
    mkdir(resdir2)
end

% Create trial events structure
if isbeh
    if isempty(protocoltag)
        if exist('TE.mat','file') ~= 0
            TE      = load('TE.mat');
        elseif exist('TE.mat','file') == 0
            behfile = dir([fullpth filesep animalID '_GamblingTask_' sessionID(1:end-1) '*.mat']);
            behfile = behfile.name;
            TE      = TE_gamblingtask(fullfile(fullpth,behfile),1);
        end
%     else % What is this for?
%         evalstr = ['TE = solo2trialevents4_pavlovian2_Sanchari_' protocoltag '([fullpth ''data_@auditory_pavlovian2_' protocoltag '_Sanchari_'' animalID2 ''_'' sessionID ''.mat'']);'];
%         eval(evalstr)
    end
    
    if isrec
        disp('You dont have recordings so far! If you do, change line 93 to your protocol');
        MakeTrialEvents2_gonogo_p(fullpth)  % synchronize
    end
end

% Update CellBase
if isrec
    addnewcells('dir', [animalID '\' sessionID '\'], 'prealign', true);
    cellids = findcell('rat',animalID,'session',sessionID);
    disp(cellids)
end

% Behavior
stg = num2str(TE.TrainingStage(1,1));
if isbeh
    switch (TE.TrainingStage(1,1))	
     case 1
        % Lick raster 
	    H = figure;
        TE.StimulusOn = zeros(1,length(TE.TrainingStage));
        viewlickRL({animalID sessionID},'TE',TE,'TriggerName','StimulusOn','SortEvent',...
            'TrialStart','eventtype','behav','ShowEvents',{{'StimulusOn'}},...
            'Partitions','#TrialType','window',[-5 5],'LickSide','all');
% 	    maximize_figure(H);
	    fnm = fullfile(resdir2,[animalID '_' sessionID '_Stg' stg '_Licks.jpg']);   % save
	    saveas(H,fnm)
	    close(H)
        
     case 2
        % Lick raster all licks, pertition trial
	    H = figure;
	    viewlickRL({animalID sessionID},'TE',TE,'TriggerName','StimulusOn','SortEvent',...
            'TrialStart','eventtype','behav','ShowEvents',{{'StimulusOn'}},...
            'Partitions','#TrialType','window',[-5 5],'LickSide','all');
% 	    maximize_figure(H);
	    fnm = fullfile(resdir2,[animalID '_' sessionID '_Stg' stg '_Licks.jpg']);   % save
	    saveas(H,fnm)
	    close(H)
      
        % Lick raster R
	    H = figure;
	    viewlickRL({animalID sessionID},'TE',TE,'TriggerName','StimulusOn','SortEvent',...
            'TrialStart','eventtype','behav','ShowEvents',{{'StimulusOn'}},...
            'Partitions','#TrialType','window',[-5 5],'LickSide','R');
% 	    maximize_figure(H);
	    fnm = fullfile(resdir2,[animalID '_' sessionID '_Stg' stg '_LicksR.jpg']);   % save
	    saveas(H,fnm)
	    close(H)
    
        % Lick raster L
	    G = figure;
	    viewlickRL({animalID sessionID},'TE',TE,'TriggerName','StimulusOn','SortEvent',...
            'TrialStart','eventtype','behav','ShowEvents',{{'StimulusOn'}},...
            'Partitions','#TrialType','window',[-5 5],'LickSide','L');
% 	    maximize_figure(G);
	    fnm = fullfile(resdir2,[animalID '_' sessionID '_Stg' stg '_LicksL.jpg']);   % save
	    saveas(G,fnm)
	    close(G)
        
        % Lick raster TrialT1
        H = figure;
	    viewlickTT({animalID sessionID},'TE',TE,'TriggerName','StimulusOn','SortEvent',...
            'TrialStart','eventtype','behav','ShowEvents',{{'StimulusOn'}},...
            'Partitions','#FirstLick','window',[-5 5],'LickSide','all','TrialType','1');
% 	    maximize_figure(H);
	    fnm = fullfile(resdir2,[animalID '_' sessionID '_Stg' stg '_TrialT1.jpg']);   % save
	    saveas(H,fnm)
	    close(H)
        
        % Lick raster TrialT2
        G = figure;
	    viewlickTT({animalID sessionID},'TE',TE,'TriggerName','StimulusOn','SortEvent',...
            'TrialStart','eventtype','behav','ShowEvents',{{'StimulusOn'}},...
            'Partitions','#FirstLick','window',[-5 5],'LickSide','all','TrialType','2');
% 	    maximize_figure(G);
	    fnm = fullfile(resdir2,[animalID '_' sessionID '_Stg' stg '_TrialT2.jpg']);   % save
	    saveas(G,fnm)
	    close(G)
        
     case 3
% 	    % Lick raster R
% 	    H = figure;
% 	    viewlickRL({animalID sessionID},'TE',TE,'TriggerName','StimulusOn','SortEvent',...
%             'TrialStart','eventtype','behav','ShowEvents',{{'StimulusOn'}},...
%             'Partitions','#TrialType','window',[-5 5],'LickSide','R');
% 	    maximize_figure(H);
% 	    fnm = fullfile(resdir2,[animalID '_' sessionID '_Stg' stg '_LicksR.jpg']);   % save
% 	    saveas(H,fnm)
% 	    close(H)
%     
%       % Lick raster L
% 	    G = figure;
% 	    viewlickRL({animalID sessionID},'TE',TE,'TriggerName','StimulusOn','SortEvent',...
%             'TrialStart','eventtype','behav','ShowEvents',{{'StimulusOn'}},...
%             'Partitions','#TrialType','window',[-5 5],'LickSide','L');
% 	    maximize_figure(G);
% 	    fnm = fullfile(resdir2,[animalID '_' sessionID '_Stg' stg '_LicksL.jpg']);   % save
% 	    saveas(G,fnm)
% 	    close(G)
         
        % Lick raster TrialT1
        H = figure;
	    viewlickTT({animalID sessionID},'TE',TE,'TriggerName','StimulusOn','SortEvent',...
            'TrialStart','eventtype','behav','ShowEvents',{{'StimulusOn'}},...
            'Partitions','#FirstLick','window',[-5 5],'LickSide','all','TrialType','1');
% 	    maximize_figure(H);
	    fnm = fullfile(resdir2,[animalID '_' sessionID '_Stg' stg '_TrialT1.jpg']);   % save
	    saveas(H,fnm)
	    close(H)
        
        % Lick raster TrialT2
        G = figure;
	    viewlickTT({animalID sessionID},'TE',TE,'TriggerName','StimulusOn','SortEvent',...
            'TrialStart','eventtype','behav','ShowEvents',{{'StimulusOn'}},...
            'Partitions','#FirstLick','window',[-5 5],'LickSide','all','TrialType','2');
% 	    maximize_figure(G);
	    fnm = fullfile(resdir2,[animalID '_' sessionID '_Stg' stg '_TrialT2.jpg']);   % save
	    saveas(G,fnm)
	    close(G)
    end   
end

% % Response profiles 
% if isbeh && isrec
%     
%     % Prealign spikes for trial events
%     problem_behav_cellid = [];
%     for iC = 1:length(cellids)
%         cellid = cellids(iC);
%         try
%             prealignSpikes(cellid,'FUNdefineEventsEpochs',@defineEventsEpochs_pavlovian,'filetype','event','ifsave',1,'ifappend',0)
%         catch
%             disp('Error in prealignSpikes.');
%             problem_behav_cellid = [problem_behav_cellid cellid];
%         end
%     end
%     
%     % Aligning to stimulus onset;
%     for k = 1:length(cellids)
%         G = figure;
%         pause(0.01)
%         viewcell2b(cellids(k),'TriggerName','StimulusOn','SortEvent','TrialStart','sigma', 0.07,'eventtype','behav','ShowEvents',{{'DeliverAllFeedback'}},'Partitions','#Reward','window',[-3 3])
%         maximize_figure(G)
%         
%         cellidt = cellids{k};
%         cellidt(cellidt=='.') = '_';
%         fnm = fullfile(resdir,[cellidt '_StimulusOn.jpg']);   % save
%         saveas(G,fnm)
%         close(G)
%     end
%     
%     for k = 1:length(cellids)
%         G = figure;
%         pause(0.01)
%         viewcell2b(cellids(k),'TriggerName','DeliverAllFeedback','SortEvent','TrialStart','sigma', 0.07,'eventtype','behav','ShowEvents',{{'StimulusOn'}},'Partitions','#Reward','window',[-3 3])
%         maximize_figure(G)
%         
%         cellidt = cellids{k};
%         cellidt(cellidt=='.') = '_';
%         fnm = fullfile(resdir,[cellidt '_Reward.jpg']);   % save
%         saveas(G,fnm)
%         close(G)
%     end
%     
%         
%     for k = 1:length(cellids)
%         G = figure;
%         pause(0.01)
%         viewcell2b(cellids(k),'TriggerName','DeliverAllFeedback','SortEvent','TrialStart','sigma', 0.07,'eventtype','behav','ShowEvents',{{'StimulusOn'}},'Partitions','#AllReward','window',[-5 5])
%         maximize_figure(G)
%         
%         cellidt = cellids{k};
%         cellidt(cellidt=='.') = '_';
%         fnm = fullfile(resdir,[cellidt '_allRPE.jpg']);   % save
%         saveas(G,fnm)
%         close(G)
%     end
%     
%     % Feedback
%     if ~isempty(nansum(TE.Hit))
%         for k = 1:length(cellids)
%             J = figure;
%             pause(0.01)
%             viewcell2b(cellids(k),'TriggerName','DeliverAllFeedback','SortEvent','TrialStart','sigma', 0.07,'eventtype','behav','ShowEvents',{{'StimulusOn'}},'Partitions','#TrialType','window',[-5 5])
%             maximize_figure(J)
%             cellidt = cellids{k};
%             cellidt(cellidt=='.') = '_';
%             fnm = fullfile(resdir,[cellidt  '_deliverallfeedback.jpg']);   % save
%             saveas(J,fnm)
%             close(J)
%         end
%         
%         if ~all(isnan(TE.Punishment))
%             for k = 1:length(cellids)
%                 J = figure;
%                 pause(0.01)
%                 viewcell2b(cellids(k),'TriggerName','DeliverAllFeedback','SortEvent','TrialStart','sigma', 0.07,'eventtype','behav','ShowEvents',{{'StimulusOn'}},'Partitions','#Punishment','window',[-5 5]);
%                 maximize_figure(J)
%                 cellidt = cellids{k};
%                 cellidt(cellidt=='.') = '_';
%                 fnm = fullfile(resdir,[cellidt '_PPE.jpg']);   % save
%                 saveas(J,fnm)
%                 close(J)
%             end
%         end
%         
%         if ~all(isnan(TE.Omission))
%             for k = 1:length(cellids)
%                 J = figure;
%                 pause(0.01)
%                 viewcell2b(cellids(k),'TriggerName','DeliverAllFeedback','SortEvent','TrialStart','sigma', 0.07,'eventtype','behav','ShowEvents',{{'StimulusOn'}},'Partitions','#Omission','window',[-5 5])
%                 maximize_figure(J)
%                 
%                 cellidt = cellids{k};
%                 cellidt(cellidt=='.') = '_';
%                 fnm = fullfile(resdir,[cellidt '_omission.jpg']);   % save
%                 saveas(J,fnm)
%                 close(J)
%             end
%         end
%         
%     end
% end
% 
% % Light effects
% if isrec && isstim
%     
%     
%     % Create stimulus events
%     MakeStimEvents2_p(fullpth,'BurstStartNttl',4)
%     
%     % Prealign spikes to stimulus events
%     problem_stim_cellid = [];
%     for iC = 1:length(cellids)
%         cellid = cellids(iC);
%         try
%             prealignSpikes(cellid,'FUNdefineEventsEpochs',@defineEventsEpochs_laserstim,'filetype','stim','ifsave',1,'ifappend',0)
%         catch
%             disp('Error in prealignSpikes.');
%             problem_stim_cellid = [problem_stim_cellid cellid];
%         end
%     end
%     
%     % View light-triggered raster and PSTH
%     TrigEvent = 'BurstOn';
%     SEvent = 'BurstOff';
%     win = [-0.2 0.5];
%     % parts = 'all';
%     parts = '#BurstNPulse';
%     dt = 0.001;
%     sigma = 0.001;
%     PSTHstd = 'on';
%     ShEvent = {{'PulseOn','PulseOff','BurstOff'}};
%     ShEvColors = hsv(length(ShEvent{1}));
%     ShEvColors = mat2cell(ShEvColors,ones(size(ShEvColors,1),1),3);
%     for iCell = 1:length(cellids)
%         cellid = cellids(iCell);
%         H = figure;
%         viewcell2b(cellid,'TriggerName',TrigEvent,'SortEvent',SEvent,'ShowEvents',ShEvent,'ShowEventsColors',{ShEvColors},...
%             'FigureNum',H,'eventtype','stim','window',win,'dt',dt,'sigma',sigma,'PSTHstd',PSTHstd,'Partitions',parts,...
%             'EventMarkerWidth',0,'PlotZeroLine','off')
%         maximize_figure(H)
%         
%         cellidt = cellid{1};
%         cellidt(cellidt=='.') = '_';
%         fnm = fullfile(resdir,[cellidt '_LS.jpg']);   % save
%         saveas(H,fnm)
%         close(H)
%     end
% end
% 
% % Cluster quality
% global CheckCluster_DrawAllWaves
% if isrec
%     CheckCluster_DrawAllWaves = true;
%     BatchSessionClust(fullpth)
%     CheckCluster_DrawAllWaves = false;
% end

end