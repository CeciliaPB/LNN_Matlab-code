function TE = trialevents_gamblingtask(filepath,ifsave)

% TRIALEVENTS_GAMBLINGTASK Creates trial events structure. The events are
% extracted from the saved behavior file. Most of the events are given
% relative to the TrialStart event.
% 
% INPUTS: 
%   - filepath: where is the data located, string. 
% Varargin: 
%   - 'ifsave': uses an addidional logical input argument do determine
%   whether to save the results. Dafault 1, to save the results.
%
% OUTPUTS: 
%   - TE: a structure containing all relevant data.
% 
% Examples: 
% TRIALEVENTS_GAMBLINGTASK(SESSPATH,IFSAVE).
%
% -------------------------------------------------------------------------
% Based on solo2trialevent_auditory_cuedoutcome.m
% By Heguedus Panna.
%
% Modified by Cecília Pardo-Bellver, 2019
% Laboratory of Network Neurophysiology
% Institute of Experimantal Medicine, Hungary.
%
% MATLAB toolbox: Signal Processing Toolbox.
% -------------------------------------------------------------------------

if nargin < 2
    ifsave = 1;
end

% Determine filepath and load data
load(filepath)

% Session ID (date and 'a' or 'b')
[folder, sessionID] = fileparts(cd);

% No. of trials (exclude trial 1 and last trial may not have been completed
ntrials = SessionData.nTrials;

% -------------------------------------------------------------------------
% Preallocation of variables
% Training session parameters
[TE.NTrials, TE.TrialStart, TE.TrialEnd, TE.TrialType] = deal(nan(1,ntrials));
TE.TrainingStage = nan(1,ntrials);
TE.sessionID     = cell(1,ntrials);
TE.Animalname    = cell(1,ntrials);

% ITI
TE.ITIBeginning = nan(1,ntrials); % Begining of ITI
TE.ITIEnd       = nan(1,ntrials);
TE.ITIDuration  = nan(1,ntrials);
TE.DeliverFeedback    = nan(1,ntrials); % Relevant outcomes: reward, punish
TE.DeliverAllFeedback = nan(1,ntrials); % All outcomes: reward, punish, omission
TE.Feedback     = nan(1,ntrials); 
TE.Omission     = nan(1,ntrials); 
TE.Punishment   = nan(1,ntrials);
TE.SafeChoice   = nan(1,ntrials);
TE.RiskyChoice  = nan(1,ntrials);

% Stimulation = sound delivery
TE.SoundFrequency1  = nan(1,ntrials);
TE.SoundFrequency2  = nan(1,ntrials);
TE.StimulusOn       = nan(1,ntrials);
TE.StimulusOff      = nan(1,ntrials);
TE.StimulusDuration = nan(1,ntrials); % stimulus duration

% Response window
TE.Delay = nan(1,ntrials); % Delay period after the stimulus delivery
TE.ResponseWindow = nan(1,ntrials); % StimulusDuration + Delay
TE.ResponseWindowEnd = nan(1,ntrials);

% Reaction
TE.LickIn   = cell(1,ntrials); % breaking the beam either side
TE.LickOut  = cell(1,ntrials); % end of breaking the beam either side
TE.LickRIn  = cell(1,ntrials); % breaking the beam of the water port R
TE.LickROut = cell(1,ntrials); % end of breaking the beam R
TE.LickLIn  = cell(1,ntrials); % breaking the beam of the water port L
TE.LickLOut = cell(1,ntrials); % end of breaking the beam L
TE.RT       = nan(1,ntrials); % Not used, just leave them
TE.RTr      = nan(1,ntrials); % 1st lick to R
TE.RTl      = nan(1,ntrials); % 1st lick to L
TE.GoRT     = nan(1,ntrials); % Not used, just leave them
TE.NoGoRT   = nan(1,ntrials); % Not used, just leave them
TE.Tup      = cell(1,ntrials); % Timeups
TE.PunishValveTime          = nan(1,ntrials); % Stage 3
TE.StartRewardSafeValveTime = nan(1,ntrials); % Stage 3
TE.EndRewardSafeValveTime   = nan(1,ntrials); % Stage 3
TE.StartRewardRiskValveTime = nan(1,ntrials); % Stage 3
TE.EndRewardRiskValveTime   = nan(1,ntrials); % Stage 3
TE.StartPunishValveTime     = nan(1,ntrials); % Stage 3 
TE.EndPunishValveTime       = nan(1,ntrials); % Stage 3
TE.StartRewardValveTime     = nan(1,ntrials); % Stage 1 & 2
TE.EndRewardValveTime       = nan(1,ntrials); % Stage 1 & 2

% Outcomes
TE.Hit              = nan(1,ntrials); % lick + reward
TE.AllReward        = nan(1,ntrials); % lick + reward, partitioned to TrialTypes
TE.Reward           = nan(1,ntrials); % lick + reward, partitioned to TrialTypes
TE.FirstLick        = nan(1,ntrials); % lick 1 = R ; 2 = L
TE.CorrectRejection = nan(1,ntrials); % no lick + punishment
TE.FalseAlarm       = nan(1,ntrials); % lick + punishment
TE.Punishment       = nan(1,ntrials); % punishment
TE.Omission         = nan(1,ntrials); % omission
TE.Miss             = nan(1,ntrials); % no lick + missed reward
TE.LickOmission     = nan(1,ntrials); % omission + lick
TE.NoLickOmission   = nan(1,ntrials); % omission + no lick

% Sync correction (state machine bug)
TE.NeedsSyncCorrection = nan(1,ntrials);

% Animal name--------------------------------------------------------------
Animalname = filepath(1:(length(filepath)-33));

% Defining parameters trialwise--------------------------------------------
for currentTrial = 1:ntrials
    TE.sessionID (1,currentTrial)       = cellstr(sessionID);
    TE.NTrials(1,currentTrial)          = ntrials; % Saving number of trials
    TE.Animalname(1,currentTrial)       = cellstr(Animalname);
    TE.TrialStart(1,currentTrial)       = SessionData.TrialStartTimestamp(currentTrial); % Start of each trial
    TE.TrainingStage(1,currentTrial)    = SessionData.TrialSettings(1,currentTrial).TrainingStage;
    TE.SoundFrequency1(1,currentTrial)  = SessionData.TrialSettings(1,currentTrial).GUI.SinWavekHz1;
    TE.SoundFrequency2(1,currentTrial)  = SessionData.TrialSettings(1,currentTrial).GUI.SinWavekHz2;

    %ITI
    if isfield(SessionData.RawEvents.Trial{1,currentTrial}.States,'ITI')
        TE.ITIBeginning(1,currentTrial) = SessionData.RawEvents.Trial{1,currentTrial}.States.ITI(1);
        TE.ITIEnd(1,currentTrial)       = SessionData.RawEvents.Trial{1,currentTrial}.States.ITI(2);
        TE.ITIDuration(1,currentTrial)  = TE.ITIEnd(1,currentTrial) - TE.ITIBeginning(1,currentTrial);
    end   
end

% Trial Type --------------------------------------------------------------
% 1 = Type 1, Reward R or Large Reward R; 2 = Type 2, Reward L or Large Reward L
TE.TrialType = SessionData.TrialTypes; 

% Outcomes, depends on hit (first calculated)------------------------------
% Hit = lick + reward; FalseAlarm = lick + punishment; LickOmission =
% omission + lick; CorrectRejection = no lick + punishment; NoLickOmission
% = omission + no lick; Miss = no lick + missed reward.
for currentTrial = 1:ntrials
    if SessionData.TrialOutcome(currentTrial) == 1 || SessionData.TrialOutcome(currentTrial) == 2
        TE.Hit(1,currentTrial) = 1; % lick + reward
    elseif SessionData.TrialOutcome(currentTrial) == 5
        if isfield(SessionData.RawEvents.Trial{1,currentTrial}.Events,'GlobalTimer_End')== 0
            TE.Hit(1,currentTrial) = 1;
        elseif isfield(SessionData.RawEvents.Trial{1,currentTrial}.Events,'GlobalTimer_End')== 1
            TE.Hit(1,currentTrial) = 0; 
        end
    elseif SessionData.TrialOutcome(currentTrial) == 0
        if isnan(SessionData.RawEvents.Trial{1,currentTrial}.States.Punish)== 0
            TE.Hit(1,currentTrial) = 1;
        elseif isnan(SessionData.RawEvents.Trial{1,currentTrial}.States.Punish)== 1
            TE.Hit(1,currentTrial) = 0; 
        end
    else
    end 
    
    if TE.Hit(currentTrial) == 1 && SessionData.TrialOutcome(currentTrial) == 0
        TE.FalseAlarm(currentTrial) = 1; % lick + punishment
    elseif TE.Hit(currentTrial) == 1 && SessionData.TrialOutcome(currentTrial) == 5
        TE.LickOmission(currentTrial) = 1; % omission + lick
    elseif TE.Hit(currentTrial) == 0 && SessionData.TrialOutcome(currentTrial) == 0 
        TE.CorrectRejection(currentTrial) = 1; % no lick + punishment
    elseif TE.Hit(currentTrial) == 0 && SessionData.TrialOutcome(currentTrial) == 5
        TE.NoLickOmission(currentTrial) = 1 ; % omission + no lick
    elseif TE.Hit(currentTrial) == 0
        TE.Miss(currentTrial) = 1; % no lick + missed reward 
    end
end

% For partitioning reinforcement ------------------------------------------
for currentTrial = 1:ntrials
    
    % Reward and choice of the animal -------------------------------------
    % TE.Reward, 1 = LickR, SafeChoice; 2 = LickR, RiskChoice; 3 = LickL,
    % SafeCoice; 4 = LickL, RiskChoice.
    if SessionData.TrialOutcome(currentTrial) ==1 && SessionData.TrialTypes(currentTrial) == 1
        TE.Reward(currentTrial) = 1;
        TE.SafeChoice(currentTrial) = 1;
    elseif SessionData.TrialOutcome(currentTrial)==2 && SessionData.TrialTypes(currentTrial) == 1
        TE.Reward(currentTrial) = 2;
        TE.RiskyChoice(currentTrial) =  1;
    elseif SessionData.TrialOutcome(currentTrial)==1 && SessionData.TrialTypes(currentTrial) == 2
        TE.Reward(currentTrial) = 3;
        TE.SafeChoice(currentTrial) = 1;
    elseif SessionData.TrialOutcome(currentTrial)==2 && SessionData.TrialTypes(currentTrial) == 2
        TE.Reward(currentTrial) = 4;
        TE.RiskyChoice(currentTrial) = 1;
    end
        
    % TE.Punishment, not sure what 1 and 2 mean ---------------------------
    if TE.FalseAlarm(currentTrial) == 1  && SessionData.TrialTypes(currentTrial) ==1
        TE.Punishment(currentTrial) = 1;
    elseif TE.CorrectRejection(currentTrial) == 1  && SessionData.TrialTypes(currentTrial) ==1
        TE.Punishment(currentTrial) = 1;
    elseif TE.FalseAlarm(currentTrial) == 1  && SessionData.TrialTypes(currentTrial) ==2
        TE.Punishment(currentTrial) = 2;
    elseif TE.CorrectRejection(currentTrial) == 1  && SessionData.TrialTypes(currentTrial) ==2
        TE.Punishment(currentTrial) = 2;
    end
    
    % TE.Omission, 1 = Omissions in TType 1; 2 = Omissions in TType 2 -----
    if TE.NoLickOmission(currentTrial) == 1 || TE.LickOmission(currentTrial) == 1 && SessionData.TrialTypes(currentTrial) == 1
        TE.Omission(currentTrial) = 1;
    elseif TE.NoLickOmission(currentTrial) == 1 || TE.LickOmission(currentTrial) == 1 && SessionData.TrialTypes(currentTrial) == 2
        TE.Omission(currentTrial) = 2;
    end
    
%     % Not needed
%     if SessionData.TrialOutcome(currentTrial) ==1 && SessionData.TrialTypes(currentTrial) ==1
%         TE.AllReward(currentTrial) = 1;
%     elseif SessionData.TrialOutcome(currentTrial)==1 && SessionData.TrialTypes(currentTrial) == 2
%         TE.AllReward(currentTrial) = 2;
%     elseif SessionData.TrialOutcome(currentTrial) ==2 && SessionData.TrialTypes(currentTrial) ==1
%         TE.AllReward(currentTrial) = 1;
%     elseif SessionData.TrialOutcome(currentTrial)==2 && SessionData.TrialTypes(currentTrial) == 2
%         TE.AllReward(currentTrial) = 2;
%     end
    
end
       
% Feedback ----------------------------------------------------------------
% Feedback: 1 = Reward; 2 = Punishment; 3 = Omission.
for currentTrial = 1:ntrials
    if SessionData.TrialOutcome(currentTrial) == 1 || SessionData.TrialOutcome(currentTrial)== 2
        TE.Feedback(currentTrial) = 1; % Reward
    elseif TE.FalseAlarm(currentTrial) == 1 || TE.CorrectRejection(currentTrial) == 1
        TE.Feedback(currentTrial) = 2; % Punishment
    elseif TE.NoLickOmission(currentTrial) == 1 || TE.LickOmission(currentTrial) == 1
        TE.Feedback(currentTrial) = 3; % Omission
    end
end

% Stimulus variables ------------------------------------------------------
for currentTrial = 1:ntrials
    if isfield(SessionData.RawEvents.Trial{1,currentTrial}.States,'StartStimulus')   
        TE.StimulusOn(1,currentTrial) = SessionData.RawEvents.Trial{1,currentTrial}.States.StartStimulus(1:1);
        TE.Delay(1,currentTrial) = SessionData.RawEvents.Trial{1, currentTrial}.States.Delay(1:1)-SessionData.RawEvents.Trial{1,currentTrial}.States.Delay(2);
    end
    
    %DeliverFeedback ------------------------------------------------------
     switch (TE.TrainingStage(1,1))	
      case 1 % TrainingStage 1
        if isfield(SessionData.RawEvents.Trial{1,currentTrial}.States,'StartStimulus')   % free water trials don't have these
            if TE.Feedback(currentTrial) == 2 % lick& no lick punishment
                TE.DeliverFeedback(currentTrial) = SessionData.RawEvents.Trial{1, currentTrial}.States.Punish(1:1) ;
            elseif TE.Feedback(currentTrial) == 1 % lick & reward
                if isnan(SessionData.RawEvents.Trial{1,currentTrial}.States.Reward) == 0
                    TE.DeliverFeedback(currentTrial) = SessionData.RawEvents.Trial{1, currentTrial}.States.Reward(1:1);
                end
            end
            if TE.Miss(currentTrial) == 1    %no lick & no reward
                TE.DeliverFeedback(currentTrial) = SessionData.RawEvents.Trial{1, currentTrial}.States.PostUS(1:1);
            end

            %DeliverAllFeedback -----------------------------------------------
            if TE.Feedback(currentTrial) == 2  % punishment
                TE.DeliverAllFeedback(currentTrial) = SessionData.RawEvents.Trial{1, currentTrial}.States.Punish(1:1) ;
            elseif TE.Feedback(currentTrial) == 1 % lick & reward
                if isnan(SessionData.RawEvents.Trial{1,currentTrial}.States.Reward) == 0
                    TE.DeliverFeedback(currentTrial) = SessionData.RawEvents.Trial{1, currentTrial}.States.Reward(1:1);
                end
            elseif TE.Feedback(currentTrial) == 3
                TE.DeliverAllFeedback(currentTrial) = SessionData.RawEvents.Trial{1, currentTrial}.States.PostUS(1:1);
            end
        end
          
      case 2 % TrainingStage 2
        if isfield(SessionData.RawEvents.Trial{1,currentTrial}.States,'StartStimulus')   % free water trials don't have these
            if TE.Feedback(currentTrial) == 2 % lick& no lick punishment
                TE.DeliverFeedback(currentTrial) = SessionData.RawEvents.Trial{1, currentTrial}.States.Punish(1:1) ;
            elseif TE.Feedback(currentTrial) == 1 % lick & reward
                if isnan(SessionData.RawEvents.Trial{1,currentTrial}.States.Reward) == 0
                    TE.DeliverFeedback(currentTrial) = SessionData.RawEvents.Trial{1, currentTrial}.States.Reward(1:1);
                end
            end
            if TE.Miss(currentTrial) == 1    %no lick & no reward
                TE.DeliverFeedback(currentTrial) = SessionData.RawEvents.Trial{1, currentTrial}.States.PostUS(1:1);
            end

            %DeliverAllFeedback -----------------------------------------------
            if TE.Feedback(currentTrial) == 2  % punishment
                TE.DeliverAllFeedback(currentTrial) = SessionData.RawEvents.Trial{1, currentTrial}.States.Punish(1:1) ;
            elseif TE.Feedback(currentTrial) == 1 % lick & reward
                if isnan(SessionData.RawEvents.Trial{1,currentTrial}.States.Reward) == 0
                    TE.DeliverFeedback(currentTrial) = SessionData.RawEvents.Trial{1, currentTrial}.States.Reward(1:1);
                end
            elseif TE.Feedback(currentTrial) == 3
                TE.DeliverAllFeedback(currentTrial) = SessionData.RawEvents.Trial{1, currentTrial}.States.PostUS(1:1);
            end
        end
        
      case 3 % TrainingStage 3       
        if isfield(SessionData.RawEvents.Trial{1,currentTrial}.States,'StartStimulus')   % free water trials don't have these
            if TE.Feedback(currentTrial) == 2 % lick& no lick punishment
                TE.DeliverFeedback(currentTrial) = SessionData.RawEvents.Trial{1, currentTrial}.States.Punish(1:1) ;
            elseif TE.Feedback(currentTrial) == 1 % lick & reward
                if isnan(SessionData.RawEvents.Trial{1,currentTrial}.States.RewardRisk) == 0
                    TE.DeliverFeedback(currentTrial) = SessionData.RawEvents.Trial{1, currentTrial}.States.RewardRisk(1:1);
                elseif isnan(SessionData.RawEvents.Trial{1,currentTrial}.States.RewardSafe) == 0
                    TE.DeliverFeedback(currentTrial) = SessionData.RawEvents.Trial{1, currentTrial}.States.RewardSafe(1:1);
                end
            end
            if TE.Miss(currentTrial) == 1    %no lick & no reward
                TE.DeliverFeedback(currentTrial) = SessionData.RawEvents.Trial{1, currentTrial}.States.PostUS(1:1);
            end

            %DeliverAllFeedback -----------------------------------------------
            if TE.Feedback(currentTrial) == 2  % punishment
                TE.DeliverAllFeedback(currentTrial) = SessionData.RawEvents.Trial{1, currentTrial}.States.Punish(1:1) ;
            elseif TE.Feedback(currentTrial) == 1 % lick & reward
                if isnan(SessionData.RawEvents.Trial{1,currentTrial}.States.RewardRisk) == 0
                    TE.DeliverFeedback(currentTrial) = SessionData.RawEvents.Trial{1, currentTrial}.States.RewardRisk(1:1);
                elseif isnan(SessionData.RawEvents.Trial{1,currentTrial}.States.RewardSafe) == 0
                    TE.DeliverFeedback(currentTrial) = SessionData.RawEvents.Trial{1, currentTrial}.States.RewardSafe(1:1);
                end
            elseif TE.Feedback(currentTrial) == 3
                TE.DeliverAllFeedback(currentTrial) = SessionData.RawEvents.Trial{1, currentTrial}.States.PostUS(1:1);
            end
        end
     end
    
end

% Reaction ----------------------------------------------------------------
for currentTrial = 1:(ntrials-1)
    if TE.Reward(currentTrial) == 1 || TE.Reward(currentTrial) == 2 % 1 = LickR, SafeChoice; 2 = LickR
        TE.LickRIn{1, currentTrial}  = SessionData.RawEvents.Trial{1, currentTrial}.Events.Port1In; % breaking the beam of the water port in headfixed protocol
        if isfield(SessionData.RawEvents.Trial{1, currentTrial}.Events, 'Port1Out')
            TE.LickROut{1, currentTrial} = SessionData.RawEvents.Trial{1, currentTrial}.Events.Port1Out;
        else
        end
        TE.LickIn{1,currentTrial}    = TE.LickRIn{1,currentTrial};
        TE.LickOut{1,currentTrial}   = TE.LickROut{1,currentTrial};
        if isfield(SessionData.RawEvents.Trial{1,currentTrial}.Events,'Port3In')
            TE.LickLIn{1, currentTrial}  = SessionData.RawEvents.Trial{1, currentTrial}.Events.Port3In;% breaking the beam of the water port in headfixed protocol
            if isfield(SessionData.RawEvents.Trial{1, currentTrial}.Events, 'Port3Out')
            TE.LickLOut{1, currentTrial} = SessionData.RawEvents.Trial{1, currentTrial}.Events.Port3Out;
            else
            end
            TE.LickIn{1,currentTrial}    = [TE.LickLIn{1,currentTrial}, TE.LickRIn{1,currentTrial}];
            TE.LickOut{1,currentTrial}   = [TE.LickLOut{1,currentTrial}, TE.LickROut{1,currentTrial}];
        end
    elseif TE.Reward(currentTrial) == 3 || TE.Reward(currentTrial) == 4 %  3 = LickL, SafeCoice; 4 = LickL, RiskChoice.
        TE.LickLIn{1, currentTrial}  = SessionData.RawEvents.Trial{1, currentTrial}.Events.Port3In;% breaking the beam of the water port in headfixed protocol
        if isfield(SessionData.RawEvents.Trial{1, currentTrial}.Events, 'Port3Out')
            TE.LickLOut{1, currentTrial} = SessionData.RawEvents.Trial{1, currentTrial}.Events.Port3Out;
        else
        end
        TE.LickIn{1,currentTrial}    = TE.LickLIn{1,currentTrial};
        TE.LickOut{1,currentTrial}   = TE.LickLOut{1,currentTrial};
        if isfield(SessionData.RawEvents.Trial{1,currentTrial}.Events,'Port1In')
            TE.LickRIn{1, currentTrial}  = SessionData.RawEvents.Trial{1, currentTrial}.Events.Port1In;% breaking the beam of the water port in headfixed protocol
            if isfield(SessionData.RawEvents.Trial{1, currentTrial}.Events, 'Port1Out')
            TE.LickROut{1, currentTrial} = SessionData.RawEvents.Trial{1, currentTrial}.Events.Port1Out;
            else
            end
            TE.LickIn{1,currentTrial}    = [TE.LickLIn{1,currentTrial}, TE.LickRIn{1,currentTrial}];
            TE.LickOut{1,currentTrial}   = [TE.LickLOut{1,currentTrial}, TE.LickROut{1,currentTrial}];
        end
    elseif TE.Feedback(currentTrial) == 3 % Omissions with lick
        if isfield(SessionData.RawEvents.Trial{1,currentTrial}.Events,'Port3In')
            TE.LickLIn{1, currentTrial}  = SessionData.RawEvents.Trial{1, currentTrial}.Events.Port3In;% breaking the beam of the water port in headfixed protocol
            if isfield(SessionData.RawEvents.Trial{1, currentTrial}.Events, 'Port3Out')
            TE.LickLOut{1, currentTrial} = SessionData.RawEvents.Trial{1, currentTrial}.Events.Port3Out;
            else
            end
            TE.LickIn{1,currentTrial}    = TE.LickLIn{1,currentTrial};
            TE.LickOut{1,currentTrial}   = TE.LickLOut{1,currentTrial};
            if isfield(SessionData.RawEvents.Trial{1,currentTrial}.Events,'Port1In')
                TE.LickRIn{1, currentTrial}  = SessionData.RawEvents.Trial{1, currentTrial}.Events.Port1In;% breaking the beam of the water port in headfixed protocol
                if isfield(SessionData.RawEvents.Trial{1, currentTrial}.Events, 'Port1Out')
                TE.LickROut{1, currentTrial} = SessionData.RawEvents.Trial{1, currentTrial}.Events.Port1Out;
                else
                end
                TE.LickIn{1,currentTrial}    = [TE.LickLIn{1,currentTrial}, TE.LickRIn{1,currentTrial}];
                TE.LickOut{1,currentTrial}   = [TE.LickLOut{1,currentTrial}, TE.LickROut{1,currentTrial}];
            end
        end
        if isfield(SessionData.RawEvents.Trial{1,currentTrial}.Events,'Port1In')
            TE.LickRIn{1, currentTrial}  = SessionData.RawEvents.Trial{1, currentTrial}.Events.Port1In;% breaking the beam of the water port in headfixed protocol
            if isfield(SessionData.RawEvents.Trial{1, currentTrial}.Events, 'Port1Out')
            TE.LickROut{1, currentTrial} = SessionData.RawEvents.Trial{1, currentTrial}.Events.Port1Out;
            else
            end
            TE.LickIn{1,currentTrial}    = TE.LickRIn{1,currentTrial};
            TE.LickOut{1,currentTrial}   = TE.LickROut{1,currentTrial};
            if isfield(SessionData.RawEvents.Trial{1,currentTrial}.Events,'Port3In')
                TE.LickLIn{1, currentTrial}  = SessionData.RawEvents.Trial{1, currentTrial}.Events.Port3In;% breaking the beam of the water port in headfixed protocol
                if isfield(SessionData.RawEvents.Trial{1, currentTrial}.Events, 'Port3Out')
                TE.LickLOut{1, currentTrial} = SessionData.RawEvents.Trial{1, currentTrial}.Events.Port3Out;
                else
                end
                TE.LickIn{1,currentTrial}    = [TE.LickLIn{1,currentTrial}, TE.LickRIn{1,currentTrial}];
                TE.LickOut{1,currentTrial}   = [TE.LickLOut{1,currentTrial}, TE.LickROut{1,currentTrial}];
            end
        end
    end
    
    if ~isempty(TE.LickRIn{1,currentTrial}) 
        rellicktimesR = TE.LickRIn{1,currentTrial}-TE.StimulusOn(1,currentTrial);
        firstlickinxR = find(rellicktimesR>0,1,'first');
        if ~isempty(firstlickinxR)   % only if there were licks after stimulus
            TE.RTr(1,currentTrial) = rellicktimesR(firstlickinxR);
        end
    end
    if ~isempty(TE.LickLIn{1,currentTrial}) 
        rellicktimesL = TE.LickLIn{1,currentTrial}-TE.StimulusOn(1,currentTrial);
        firstlickinxL = find(rellicktimesL>0,1,'first');
        if ~isempty(firstlickinxL)   % only if there were licks after stimulus
        TE.RTl(1,currentTrial) = rellicktimesL(firstlickinxL);
        end
    end
    
    if TE.RTr(1,currentTrial) > TE.RTl(1,currentTrial)
        TE.RT = TE.RTl;
    elseif TE.RTr(1,currentTrial) < TE.RTl(1,currentTrial)
        TE.RT = TE.RTr;
    else
    end
    
    
end
[TE.GoRT, TE.NoGoRT]        = deal(TE.RT);
TE.GoRT(TE.TrialType~=1)    = NaN;   % 'reaction time' in the likely rewarded trials
TE.NoGoRT(TE.TrialType~=2)  = NaN;   % 'reaction time' in the likely punished trials

for currentTrial = 1:ntrials
    switch (TE.TrainingStage(1,1))	
      case 1 % TrainingStage 1
        TE.Tup{1, currentTrial} = SessionData.RawEvents.Trial{1, currentTrial}.Events.Tup;
        TE.StartRewardValveTime(1,currentTrial) = SessionData.RawEvents.Trial{1, currentTrial}.States.Reward(1:1);
        TE.EndRewardValveTime(1,currentTrial)   = SessionData.RawEvents.Trial{1, currentTrial}.States.Reward(2);
          
      case 2 % TrainingStage 2
        TE.Tup{1, currentTrial} = SessionData.RawEvents.Trial{1, currentTrial}.Events.Tup;
        TE.StartRewardValveTime(1,currentTrial) = SessionData.RawEvents.Trial{1, currentTrial}.States.Reward(1:1);
        TE.EndRewardValveTime(1,currentTrial)   = SessionData.RawEvents.Trial{1, currentTrial}.States.Reward(2);
          
      case 3 % TrainingStage 3
        TE.Tup{1, currentTrial} = SessionData.RawEvents.Trial{1, currentTrial}.Events.Tup;
        TE.StartRewardSafeValveTime(1,currentTrial) = SessionData.RawEvents.Trial{1, currentTrial}.States.RewardSafe(1:1);
        TE.EndRewardSafeValveTime(1,currentTrial)   = SessionData.RawEvents.Trial{1, currentTrial}.States.RewardSafe(2);
        TE.StartRewardRiskValveTime(1,currentTrial) = SessionData.RawEvents.Trial{1, currentTrial}.States.RewardRisk(1:1);
        TE.EndRewardRiskValveTime(1,currentTrial)   = SessionData.RawEvents.Trial{1, currentTrial}.States.RewardRisk(2);
        TE.StartPunishValveTime(1,currentTrial)     = SessionData.RawEvents.Trial{1, currentTrial}.States.Punish(1:1);
        TE.EndPunishValveTime(1,currentTrial)       = SessionData.RawEvents.Trial{1, currentTrial}.States.Punish(2);
    end
end

% Sync correction
for currentTrial = 1:ntrials
    if isfield(SessionData.RawEvents.Trial{1,currentTrial}.States,'DeliverStimulus')   % TTL was sent with 25 ms delay by mistake
        TE.NeedsSyncCorrection(currentTrial) = 1;
    end
end

% Save
if ifsave == 1
    
    % For each session, save the individual variables from the extracted TE
    if ~exist([folder filesep sessionID],'dir') == 1
        mkdir([folder filesep sessionID]);
    end
    TEdir = [folder filesep sessionID];
    savename = 'TE.mat';
    save([TEdir filesep savename],'-struct','TE')
    
    % this save process saves individual variables in the file; can put into a
    % struct format by writing TE=load(TE_hr012_100718a.mat), and then can call
    % up TE.sessionID etc.
end

disp('TE analysis complete');
end

% -------------------------------------------------------------------------
function TE2 = shortenTE(TE2,shinx)

% Eliminate behavioral trials
fnm = fieldnames(TE2);
for k = 1:length(fieldnames(TE2))
    TE2.(fnm{k}) = TE2.(fnm{k})(shinx);
end
end