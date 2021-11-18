function TE_behaviour = TE_gamblingtask(filepath,ifsave)

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
%   - TE_behaviour: a structure containing all relevant data.
% 
% Examples: 
% TE_GAMBLINGTASK(SESSPATH,IFSAVE).
%
% -------------------------------------------------------------------------
% Based on solo2trialevent_auditory_cuedoutcome.m by Heguedus Panna.
%
% Cecília Pardo-Bellver, 2019
% Laboratory of Network Neurophysiology
% Institute of Experimental Medicine, Hungary.
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

% Preallocation of variables-----------------------------------------------
% Training session parameters
[TE_behaviour.NTrials, TE_behaviour.TrialStart, TE_behaviour.TrialEnd, TE_behaviour.TrialType] = deal(nan(1,ntrials));
TE_behaviour.TrainingStage = nan(1,ntrials);
TE_behaviour.sessionID     = cell(1,ntrials);
TE_behaviour.Animalname    = cell(1,ntrials);

% ITI
TE_behaviour.ITIBeginning = nan(1,ntrials); % Begining of ITI
TE_behaviour.ITIEnd       = nan(1,ntrials);
TE_behaviour.ITIDuration  = nan(1,ntrials);

% Stimulation, sound delivery
TE_behaviour.SoundFrequency1  = nan(1,ntrials);
TE_behaviour.SoundFrequency2  = nan(1,ntrials);
TE_behaviour.StimulusOn       = nan(1,ntrials);
TE_behaviour.StimulusOff      = nan(1,ntrials);
TE_behaviour.StimulusDuration = nan(1,ntrials); % stimulus duration

% Outcomes
TE_behaviour.Hit              = nan(1,ntrials); % lick + reward
TE_behaviour.Reward           = nan(1,ntrials); % lick + reward, partitioned to TrialTypes
TE_behaviour.CorrectRejection = nan(1,ntrials); % no lick + punishment
TE_behaviour.FalseAlarm       = nan(1,ntrials); % lick + punishment
TE_behaviour.Punishment       = nan(1,ntrials); % punishment
TE_behaviour.Omission         = nan(1,ntrials); % omission
TE_behaviour.Miss             = nan(1,ntrials); % no lick + missed reward
TE_behaviour.LickOmission     = nan(1,ntrials); % omission + lick
TE_behaviour.NoLickOmission   = nan(1,ntrials); % omission + no lick
TE_behaviour.SafeChoice       = nan(1,ntrials); % Safe choice
TE_behaviour.RiskyChoice      = nan(1,ntrials); % Risk choice

% Reaction & Stimulus depending variables
TE_behaviour.DeliverFeedback    = nan(1,ntrials); % Time of reward, punish
TE_behaviour.DeliverAllFeedback = nan(1,ntrials); % Time of reward, punish, omission
TE_behaviour.Feedback  = nan(1,ntrials);  % 1 = Reward; 2 = Punishment; 3 = Omission.
TE_behaviour.LickIn    = cell(1,ntrials); % breaking the beam either side
TE_behaviour.LickOut   = cell(1,ntrials); % end of breaking the beam either side
TE_behaviour.LickRIn   = cell(1,ntrials); % breaking the beam of the water port R
TE_behaviour.LickROut  = cell(1,ntrials); % end of breaking the beam R
TE_behaviour.LickLIn   = cell(1,ntrials); % breaking the beam of the water port L
TE_behaviour.LickLOut  = cell(1,ntrials); % end of breaking the beam L
TE_behaviour.FirstLick = nan(1,ntrials); % lick 1 = R ; 2 = L
TE_behaviour.RT        = nan(1,ntrials); % Timing first lick
TE_behaviour.RTr       = nan(1,ntrials); % 1st lick to R
TE_behaviour.RTl       = nan(1,ntrials); % 1st lick to L
TE_behaviour.GoRT      = nan(1,ntrials); % ! Not used, just leave them
TE_behaviour.NoGoRT    = nan(1,ntrials); % ! Not used, just leave them
TE_behaviour.Tup       = cell(1,ntrials); % Timeups
TE_behaviour.PunishValveTime          = nan(1,ntrials); % Stage 3
TE_behaviour.StartRewardSafeValveTime = nan(1,ntrials); % Stage 3
TE_behaviour.EndRewardSafeValveTime   = nan(1,ntrials); % Stage 3
TE_behaviour.StartRewardRiskValveTime = nan(1,ntrials); % Stage 3
TE_behaviour.EndRewardRiskValveTime   = nan(1,ntrials); % Stage 3
TE_behaviour.StartPunishValveTime     = nan(1,ntrials); % Stage 3 
TE_behaviour.EndPunishValveTime       = nan(1,ntrials); % Stage 3
TE_behaviour.StartRewardValveTime     = nan(1,ntrials); % Stage 1 & 2
TE_behaviour.EndRewardValveTime       = nan(1,ntrials); % Stage 1 & 2

% Response window
TE_behaviour.Delay = nan(1,ntrials); % Delay period after the stimulus delivery
TE_behaviour.ResponseWindow    = nan(1,ntrials); % StimulusDuration + Delay
TE_behaviour.ResponseWindowEnd = nan(1,ntrials);

% Sync correction (state machine bug) -------------------------------------
TE_behaviour.NeedsSyncCorrection = nan(1,ntrials);

% Animal name--------------------------------------------------------------
Animalname = filepath(1:(length(filepath)-33));

% Defining parameters trialwise--------------------------------------------
for currentTrial = 1:ntrials
    
    TE_behaviour.sessionID (1,currentTrial)       = cellstr(sessionID);
    TE_behaviour.NTrials(1,currentTrial)          = ntrials; % Saving number of trials
    TE_behaviour.Animalname(1,currentTrial)       = cellstr(Animalname);
    TE_behaviour.TrialStart(1,currentTrial)       = SessionData.TrialStartTimestamp(currentTrial); % Start of each trial
    TE_behaviour.TrainingStage(1,currentTrial)    = SessionData.TrialSettings(1,currentTrial).TrainingStage;
    TE_behaviour.SoundFrequency1(1,currentTrial)  = SessionData.TrialSettings(1,currentTrial).GUI.SinWavekHz1;
    TE_behaviour.SoundFrequency2(1,currentTrial)  = SessionData.TrialSettings(1,currentTrial).GUI.SinWavekHz2;

    % ITI 
    if isfield(SessionData.RawEvents.Trial{1,currentTrial}.States,'ITI')
        TE_behaviour.ITIBeginning(1,currentTrial) = SessionData.RawEvents.Trial{1,currentTrial}.States.ITI(1);
        TE_behaviour.ITIEnd(1,currentTrial)       = SessionData.RawEvents.Trial{1,currentTrial}.States.ITI(2);
        TE_behaviour.ITIDuration(1,currentTrial)  = TE_behaviour.ITIEnd(1,currentTrial) - TE_behaviour.ITIBeginning(1,currentTrial);
    end   
end

% Trial Type --------------------------------------------------------------
% 1 = Type 1, Reward R or Large Reward R; 2 = Type 2, Reward L or Large Reward L
TE_behaviour.TrialType = SessionData.TrialTypes; 

% Outcomes ----------------------------------------------------------------
for currentTrial = 1:ntrials
    
    % Hit, FalseAlarm, LickOmission, CorrectRejection, NoLickOmission, Miss
    if      SessionData.TrialOutcome(currentTrial) == 1 || SessionData.TrialOutcome(currentTrial) == 2 
        TE_behaviour.Hit(1,currentTrial) = 1; % lick + reward
    elseif  SessionData.TrialOutcome(currentTrial) == 3
        TE_behaviour.FalseAlarm(currentTrial) = 1; % lick + punishment
    elseif  SessionData.TrialOutcome(currentTrial) == 4
        TE_behaviour.LickOmission(currentTrial) = 1; % omission + lick
    elseif  SessionData.TrialOutcome(currentTrial) == 5 
        TE_behaviour.CorrectRejection(currentTrial) = 1; % no lick + punishment
    elseif  SessionData.TrialOutcome(currentTrial) == 6
        TE_behaviour.NoLickOmission(currentTrial) = 1 ; % omission + no lick
    elseif  SessionData.TrialOutcome(currentTrial) == 0
        TE_behaviour.Miss(currentTrial) = 1; % no lick + missed reward 
    end
   
    % Reward and choice of the animal 
    % TE.Reward, 1 = LickR, SafeChoice; 2 = LickR, RiskChoice; 3 = LickL,
    % SafeCoice; 4 = LickL, RiskChoice.
    if     SessionData.TrialOutcome(currentTrial) == 1 && SessionData.TrialTypes(currentTrial) == 1
        TE_behaviour.Reward(currentTrial) = 1;
        TE_behaviour.SafeChoice(currentTrial) = 1;
    elseif SessionData.TrialOutcome(currentTrial) == 2 && SessionData.TrialTypes(currentTrial) == 1
        TE_behaviour.Reward(currentTrial) = 2;
        TE_behaviour.RiskyChoice(currentTrial) =  1;
    elseif SessionData.TrialOutcome(currentTrial) == 1 && SessionData.TrialTypes(currentTrial) == 2
        TE_behaviour.Reward(currentTrial) = 3;
        TE_behaviour.SafeChoice(currentTrial) = 1;
    elseif SessionData.TrialOutcome(currentTrial) == 2 && SessionData.TrialTypes(currentTrial) == 2
        TE_behaviour.Reward(currentTrial) = 4;
        TE_behaviour.RiskyChoice(currentTrial) = 1;
    end
        
    % TE.Punishment, 1 for TType 1 and 2 for TType 2 
    if     TE_behaviour.FalseAlarm(currentTrial) == 1  && SessionData.TrialTypes(currentTrial) == 1
        TE_behaviour.Punishment(currentTrial) = 1;
    elseif TE_behaviour.CorrectRejection(currentTrial) == 1  && SessionData.TrialTypes(currentTrial) == 1
        TE_behaviour.Punishment(currentTrial) = 1;
    elseif TE_behaviour.FalseAlarm(currentTrial) == 1  && SessionData.TrialTypes(currentTrial) == 2 
        TE_behaviour.Punishment(currentTrial) = 2;
    elseif TE_behaviour.CorrectRejection(currentTrial) == 1  && SessionData.TrialTypes(currentTrial) == 2
        TE_behaviour.Punishment(currentTrial) = 2;
    end
    
    % TE.Omission, 1 = Omissions in TType 1; 2 = Omissions in TType 2 
    if     (TE_behaviour.NoLickOmission(currentTrial) == 1 || TE_behaviour.LickOmission(currentTrial) == 1) && SessionData.TrialTypes(currentTrial) == 1
        TE_behaviour.Omission(currentTrial) = 1;
    elseif (TE_behaviour.NoLickOmission(currentTrial) == 1 || TE_behaviour.LickOmission(currentTrial) == 1) && SessionData.TrialTypes(currentTrial) == 2
        TE_behaviour.Omission(currentTrial) = 2;
    end 
end

% Reaction & Stimulus depending variables ---------------------------------
for currentTrial = 1:ntrials
    % Feedback: 1 = Reward; 2 = Punishment; 3 = Omission.
    if     TE_behaviour.Hit(currentTrial) == 1
        TE_behaviour.Feedback(currentTrial) = 1; % Reward
    elseif TE_behaviour.FalseAlarm(currentTrial) == 1 || TE_behaviour.CorrectRejection(currentTrial) == 1
        TE_behaviour.Feedback(currentTrial) = 2; % Punishment
    elseif TE_behaviour.NoLickOmission(currentTrial) == 1 || TE_behaviour.LickOmission(currentTrial) == 1
        TE_behaviour.Feedback(currentTrial) = 3; % Omission
    end
    
    % Stimulis start, end & Delay 
    if isfield(SessionData.RawEvents.Trial{1,currentTrial}.States,'StartStimulus')   
        TE_behaviour.StimulusOn(1,currentTrial)  = SessionData.RawEvents.Trial{1,currentTrial}.States.StartStimulus(1:1);
        TE_behaviour.StimulusOff(1,currentTrial) = SessionData.RawEvents.Trial{1,currentTrial}.States.StartStimulus(1:1) + SessionData.TrialSettings(1,currentTrial).SoundDuration;
        TE_behaviour.Delay(1,currentTrial) = SessionData.RawEvents.Trial{1, currentTrial}.States.Delay(1:1)-SessionData.RawEvents.Trial{1,currentTrial}.States.Delay(2);
    end
    
    %DeliverFeedback & DeliverAllFeedback
     switch (TE_behaviour.TrainingStage(1,1))	
         
      case 1 % TrainingStage 1  % Feedback: 1 = Reward; 2 = Punishment; 3 = Omission.
        if isfield(SessionData.RawEvents.Trial{1,currentTrial}.States,'StartStimulus')   % free water trials don't have these
            if     TE_behaviour.Feedback(currentTrial) == 2 % Punishment
                TE_behaviour.DeliverFeedback(currentTrial) = SessionData.RawEvents.Trial{1, currentTrial}.States.Punish(1:1) ;
            elseif TE_behaviour.Feedback(currentTrial) == 1 % Reward
                if isnan(SessionData.RawEvents.Trial{1,currentTrial}.States.Reward) == 0
                    TE_behaviour.DeliverFeedback(currentTrial) = SessionData.RawEvents.Trial{1, currentTrial}.States.Reward(1:1);
                end
            end
            if TE_behaviour.Miss(currentTrial) == 1    % No lick & No reward
                TE_behaviour.DeliverFeedback(currentTrial) = SessionData.RawEvents.Trial{1, currentTrial}.States.PostUS(1:1);
            end

            if     TE_behaviour.Feedback(currentTrial) == 2  % Punishment
                TE_behaviour.DeliverAllFeedback(currentTrial) = SessionData.RawEvents.Trial{1, currentTrial}.States.Punish(1:1) ;
            elseif TE_behaviour.Feedback(currentTrial) == 1 % Reward
                if isnan(SessionData.RawEvents.Trial{1,currentTrial}.States.Reward) == 0
                    TE_behaviour.DeliverFeedback(currentTrial) = SessionData.RawEvents.Trial{1, currentTrial}.States.Reward(1:1);
                end
            elseif TE_behaviour.Feedback(currentTrial) == 3 % Omission
                TE_behaviour.DeliverAllFeedback(currentTrial) = SessionData.RawEvents.Trial{1, currentTrial}.States.PostUS(1:1);
            end
        end
          
      case 2 % TrainingStage 2
        if isfield(SessionData.RawEvents.Trial{1,currentTrial}.States,'StartStimulus')   % free water trials don't have these
            if     TE_behaviour.Feedback(currentTrial) == 2 % Punishment
                TE_behaviour.DeliverFeedback(currentTrial) = SessionData.RawEvents.Trial{1, currentTrial}.States.Punish(1:1) ;
            elseif TE_behaviour.Feedback(currentTrial) == 1 % Reward
                if isnan(SessionData.RawEvents.Trial{1,currentTrial}.States.Reward) == 0
                    TE_behaviour.DeliverFeedback(currentTrial) = SessionData.RawEvents.Trial{1, currentTrial}.States.Reward(1:1);
                end
            end
            if TE_behaviour.Miss(currentTrial) == 1    % No lick & No reward
                TE_behaviour.DeliverFeedback(currentTrial) = SessionData.RawEvents.Trial{1, currentTrial}.States.PostUS(1:1);
            end

            if     TE_behaviour.Feedback(currentTrial) == 2  % Punishment
                TE_behaviour.DeliverAllFeedback(currentTrial) = SessionData.RawEvents.Trial{1, currentTrial}.States.Punish(1:1) ;
            elseif TE_behaviour.Feedback(currentTrial) == 1 % Reward
                if isnan(SessionData.RawEvents.Trial{1,currentTrial}.States.Reward) == 0
                    TE_behaviour.DeliverFeedback(currentTrial) = SessionData.RawEvents.Trial{1, currentTrial}.States.Reward(1:1);
                end
            elseif TE_behaviour.Feedback(currentTrial) == 3 % Omissions
                TE_behaviour.DeliverAllFeedback(currentTrial) = SessionData.RawEvents.Trial{1, currentTrial}.States.PostUS(1:1);
            end
        end
        
      case 3 % TrainingStage 3       
        if isfield(SessionData.RawEvents.Trial{1,currentTrial}.States,'StartStimulus')   % free water trials don't have these
            if     TE_behaviour.Feedback(currentTrial) == 2 % Punishment
                TE_behaviour.DeliverFeedback(currentTrial) = SessionData.RawEvents.Trial{1, currentTrial}.States.Punish(1:1) ;
            elseif TE_behaviour.Feedback(currentTrial) == 1 % Reward
                if     isnan(SessionData.RawEvents.Trial{1,currentTrial}.States.RewardRisk) == 0
                    TE_behaviour.DeliverFeedback(currentTrial) = SessionData.RawEvents.Trial{1, currentTrial}.States.RewardRisk(1:1);
                elseif isnan(SessionData.RawEvents.Trial{1,currentTrial}.States.RewardSafe) == 0
                    TE_behaviour.DeliverFeedback(currentTrial) = SessionData.RawEvents.Trial{1, currentTrial}.States.RewardSafe(1:1);
                end
            end
            if TE_behaviour.Miss(currentTrial) == 1    % No lick & No reward
                TE_behaviour.DeliverFeedback(currentTrial) = SessionData.RawEvents.Trial{1, currentTrial}.States.PostUS(1:1);
            end

            if     TE_behaviour.Feedback(currentTrial) == 2  % Punishment
                TE_behaviour.DeliverAllFeedback(currentTrial) = SessionData.RawEvents.Trial{1, currentTrial}.States.Punish(1:1) ;
            elseif TE_behaviour.Feedback(currentTrial) == 1 % Reward
                if     isnan(SessionData.RawEvents.Trial{1,currentTrial}.States.RewardRisk) == 0
                    TE_behaviour.DeliverAllFeedback(currentTrial) = SessionData.RawEvents.Trial{1, currentTrial}.States.RewardRisk(1:1);
                elseif isnan(SessionData.RawEvents.Trial{1,currentTrial}.States.RewardSafe) == 0
                    TE_behaviour.DeliverAllFeedback(currentTrial) = SessionData.RawEvents.Trial{1, currentTrial}.States.RewardSafe(1:1);
                end
            elseif TE_behaviour.Feedback(currentTrial) == 3 % Omissions
                TE_behaviour.DeliverAllFeedback(currentTrial) = SessionData.RawEvents.Trial{1, currentTrial}.States.PostUS(1:1);
            end
        end
     end
    
end

% More reactions, last trial not included ---------------------------------
for currentTrial = 1:(ntrials-1)
    switch (TE_behaviour.TrainingStage(1,1))	
      case {1, 2} % TrainingStage1
            if     TE_behaviour.Reward(currentTrial) == 1 || TE_behaviour.Reward(currentTrial) == 2 % 1 = LickR, SafeChoice; 2 = LickR, RisckCoice.
                % R licks
                TE_behaviour.LickRIn{1, currentTrial}  = SessionData.RawEvents.Trial{1, currentTrial}.Events.Port1In; 
                if isfield(SessionData.RawEvents.Trial{1, currentTrial}.Events, 'Port1Out')
                    TE_behaviour.LickROut{1, currentTrial} = SessionData.RawEvents.Trial{1, currentTrial}.Events.Port1Out;
                else
                end
                % L licks
                try TE_behaviour.LickLIn{1, currentTrial}  = SessionData.RawEvents.Trial{1, currentTrial}.Events.Port3In;
                    if isfield(SessionData.RawEvents.Trial{1, currentTrial}.Events, 'Port3Out')
                        TE_behaviour.LickLOut{1, currentTrial} = SessionData.RawEvents.Trial{1, currentTrial}.Events.Port3Out;
                    else
                    end
                catch
                end
                TE_behaviour.LickIn{1,currentTrial}    = [TE_behaviour.LickLIn{1,currentTrial}, TE_behaviour.LickRIn{1,currentTrial}];
                TE_behaviour.LickOut{1,currentTrial}   = [TE_behaviour.LickLOut{1,currentTrial}, TE_behaviour.LickROut{1,currentTrial}];
            elseif TE_behaviour.Reward(currentTrial) == 3 || TE_behaviour.Reward(currentTrial) == 4 %  3 = LickL, SafeCoice; 4 = LickL, RiskChoice.
                % L licks
                TE_behaviour.LickLIn{1, currentTrial}  = SessionData.RawEvents.Trial{1, currentTrial}.Events.Port3In;
                if isfield(SessionData.RawEvents.Trial{1, currentTrial}.Events, 'Port3Out')
                    TE_behaviour.LickLOut{1, currentTrial} = SessionData.RawEvents.Trial{1, currentTrial}.Events.Port3Out;
                else
                end
                % R licks
                try TE_behaviour.LickRIn{1, currentTrial}  = SessionData.RawEvents.Trial{1, currentTrial}.Events.Port1In;
                    if isfield(SessionData.RawEvents.Trial{1, currentTrial}.Events, 'Port1Out')
                    TE_behaviour.LickROut{1, currentTrial} = SessionData.RawEvents.Trial{1, currentTrial}.Events.Port1Out;
                    else
                    end            
                catch
                end
                TE_behaviour.LickIn{1,currentTrial}    = [TE_behaviour.LickLIn{1,currentTrial}, TE_behaviour.LickRIn{1,currentTrial}];
                TE_behaviour.LickOut{1,currentTrial}   = [TE_behaviour.LickLOut{1,currentTrial}, TE_behaviour.LickROut{1,currentTrial}];
            elseif TE_behaviour.Feedback(currentTrial) == 3 % Omissions with lick
                if isfield(SessionData.RawEvents.Trial{1,currentTrial}.Events,'Port3In') || isfield(SessionData.RawEvents.Trial{1,currentTrial}.Events,'Port1In')
                    % L licks
                    try TE_behaviour.LickLIn{1, currentTrial}  = SessionData.RawEvents.Trial{1, currentTrial}.Events.Port3In;
                        if isfield(SessionData.RawEvents.Trial{1, currentTrial}.Events, 'Port3Out')
                        TE_behaviour.LickLOut{1, currentTrial} = SessionData.RawEvents.Trial{1, currentTrial}.Events.Port3Out;
                        else
                        end
                    catch
                    end
                    % R licks
                    try TE_behaviour.LickRIn{1, currentTrial}  = SessionData.RawEvents.Trial{1, currentTrial}.Events.Port1In;
                        if isfield(SessionData.RawEvents.Trial{1, currentTrial}.Events, 'Port1Out')
                        TE_behaviour.LickROut{1, currentTrial} = SessionData.RawEvents.Trial{1, currentTrial}.Events.Port1Out;
                        else
                        end
                    catch
                    end
                    TE_behaviour.LickIn{1,currentTrial}    = [TE_behaviour.LickLIn{1,currentTrial}, TE_behaviour.LickRIn{1,currentTrial}];
                    TE_behaviour.LickOut{1,currentTrial}   = [TE_behaviour.LickLOut{1,currentTrial}, TE_behaviour.LickROut{1,currentTrial}];
                end
            end
            
        case 3
             if     TE_behaviour.Reward(currentTrial) == 1 % 1 = LickR, SafeChoice;
                 % L licks
                TE_behaviour.LickLIn{1, currentTrial}  = SessionData.RawEvents.Trial{1, currentTrial}.Events.Port3In;
                if isfield(SessionData.RawEvents.Trial{1, currentTrial}.Events, 'Port3Out')
                    TE_behaviour.LickLOut{1, currentTrial} = SessionData.RawEvents.Trial{1, currentTrial}.Events.Port3Out;
                else
                end
                % R licks
                try TE_behaviour.LickRIn{1, currentTrial}  = SessionData.RawEvents.Trial{1, currentTrial}.Events.Port1In;
                    if isfield(SessionData.RawEvents.Trial{1, currentTrial}.Events, 'Port1Out')
                    TE_behaviour.LickROut{1, currentTrial} = SessionData.RawEvents.Trial{1, currentTrial}.Events.Port1Out;
                    else
                    end            
                catch
                end
                TE_behaviour.LickIn{1,currentTrial}    = [TE_behaviour.LickLIn{1,currentTrial}, TE_behaviour.LickRIn{1,currentTrial}];
                TE_behaviour.LickOut{1,currentTrial}   = [TE_behaviour.LickLOut{1,currentTrial}, TE_behaviour.LickROut{1,currentTrial}];
             elseif TE_behaviour.Reward(currentTrial) == 2 % 2 = LickR, RisckCoice.
                % R licks
                TE_behaviour.LickRIn{1, currentTrial}  = SessionData.RawEvents.Trial{1, currentTrial}.Events.Port1In; 
                if isfield(SessionData.RawEvents.Trial{1, currentTrial}.Events, 'Port1Out')
                    TE_behaviour.LickROut{1, currentTrial} = SessionData.RawEvents.Trial{1, currentTrial}.Events.Port1Out;
                else
                end
                % L licks
                try TE_behaviour.LickLIn{1, currentTrial}  = SessionData.RawEvents.Trial{1, currentTrial}.Events.Port3In;
                    if isfield(SessionData.RawEvents.Trial{1, currentTrial}.Events, 'Port3Out')
                        TE_behaviour.LickLOut{1, currentTrial} = SessionData.RawEvents.Trial{1, currentTrial}.Events.Port3Out;
                    else
                    end
                catch
                end
                TE_behaviour.LickIn{1,currentTrial}    = [TE_behaviour.LickLIn{1,currentTrial}, TE_behaviour.LickRIn{1,currentTrial}];
                TE_behaviour.LickOut{1,currentTrial}   = [TE_behaviour.LickLOut{1,currentTrial}, TE_behaviour.LickROut{1,currentTrial}];
             elseif TE_behaviour.Reward(currentTrial) == 3 % 3 = LickL, SafeCoice;
                 % R licks
                TE_behaviour.LickRIn{1, currentTrial}  = SessionData.RawEvents.Trial{1, currentTrial}.Events.Port1In; 
                if isfield(SessionData.RawEvents.Trial{1, currentTrial}.Events, 'Port1Out')
                    TE_behaviour.LickROut{1, currentTrial} = SessionData.RawEvents.Trial{1, currentTrial}.Events.Port1Out;
                else
                end
                % L licks
                try TE_behaviour.LickLIn{1, currentTrial}  = SessionData.RawEvents.Trial{1, currentTrial}.Events.Port3In;
                    if isfield(SessionData.RawEvents.Trial{1, currentTrial}.Events, 'Port3Out')
                        TE_behaviour.LickLOut{1, currentTrial} = SessionData.RawEvents.Trial{1, currentTrial}.Events.Port3Out;
                    else
                    end
                catch
                end
                TE_behaviour.LickIn{1,currentTrial}    = [TE_behaviour.LickLIn{1,currentTrial}, TE_behaviour.LickRIn{1,currentTrial}];
                TE_behaviour.LickOut{1,currentTrial}   = [TE_behaviour.LickLOut{1,currentTrial}, TE_behaviour.LickROut{1,currentTrial}];
             elseif TE_behaviour.Reward(currentTrial) == 4 % 4 = LickL, RiskChoice.
                % L licks
                TE_behaviour.LickLIn{1, currentTrial}  = SessionData.RawEvents.Trial{1, currentTrial}.Events.Port3In;
                if isfield(SessionData.RawEvents.Trial{1, currentTrial}.Events, 'Port3Out')
                    TE_behaviour.LickLOut{1, currentTrial} = SessionData.RawEvents.Trial{1, currentTrial}.Events.Port3Out;
                else
                end
                % R licks
                try TE_behaviour.LickRIn{1, currentTrial}  = SessionData.RawEvents.Trial{1, currentTrial}.Events.Port1In;
                    if isfield(SessionData.RawEvents.Trial{1, currentTrial}.Events, 'Port1Out')
                    TE_behaviour.LickROut{1, currentTrial} = SessionData.RawEvents.Trial{1, currentTrial}.Events.Port1Out;
                    else
                    end            
                catch
                end
                TE_behaviour.LickIn{1,currentTrial}    = [TE_behaviour.LickLIn{1,currentTrial}, TE_behaviour.LickRIn{1,currentTrial}];
                TE_behaviour.LickOut{1,currentTrial}   = [TE_behaviour.LickLOut{1,currentTrial}, TE_behaviour.LickROut{1,currentTrial}];
            elseif TE_behaviour.Feedback(currentTrial) == 3 % Omissions with lick
                if isfield(SessionData.RawEvents.Trial{1,currentTrial}.Events,'Port3In') || isfield(SessionData.RawEvents.Trial{1,currentTrial}.Events,'Port1In')
                    % L licks
                    try TE_behaviour.LickLIn{1, currentTrial}  = SessionData.RawEvents.Trial{1, currentTrial}.Events.Port3In;
                        if isfield(SessionData.RawEvents.Trial{1, currentTrial}.Events, 'Port3Out')
                        TE_behaviour.LickLOut{1, currentTrial} = SessionData.RawEvents.Trial{1, currentTrial}.Events.Port3Out;
                        else
                        end
                    catch
                    end
                    % R licks
                    try TE_behaviour.LickRIn{1, currentTrial}  = SessionData.RawEvents.Trial{1, currentTrial}.Events.Port1In;
                        if isfield(SessionData.RawEvents.Trial{1, currentTrial}.Events, 'Port1Out')
                        TE_behaviour.LickROut{1, currentTrial} = SessionData.RawEvents.Trial{1, currentTrial}.Events.Port1Out;
                        else
                        end
                    catch
                    end
                    TE_behaviour.LickIn{1,currentTrial}    = [TE_behaviour.LickLIn{1,currentTrial}, TE_behaviour.LickRIn{1,currentTrial}];
                    TE_behaviour.LickOut{1,currentTrial}   = [TE_behaviour.LickLOut{1,currentTrial}, TE_behaviour.LickROut{1,currentTrial}];
                end
            end
    end
    
    % First lick to R 
    if ~isempty(TE_behaviour.LickRIn{1,currentTrial}) 
        rellicktimesR = TE_behaviour.LickRIn{1,currentTrial}-TE_behaviour.StimulusOn(1,currentTrial);
        firstlickinxR = find(rellicktimesR>0,1,'first');
        if ~isempty(firstlickinxR)   % only if there were licks after stimulus
            TE_behaviour.RTr(1,currentTrial) = rellicktimesR(firstlickinxR);
        end
    end
    
    % First lick to L
    if ~isempty(TE_behaviour.LickLIn{1,currentTrial}) 
        rellicktimesL = TE_behaviour.LickLIn{1,currentTrial}-TE_behaviour.StimulusOn(1,currentTrial);
        firstlickinxL = find(rellicktimesL>0,1,'first');
        if ~isempty(firstlickinxL)   % only if there were licks after stimulus
        TE_behaviour.RTl(1,currentTrial) = rellicktimesL(firstlickinxL);
        end
    end
    
    % First lick: 1 == R; 2 == L
    if TE_behaviour.RTr(1,currentTrial) > TE_behaviour.RTl(1,currentTrial)
        TE_behaviour.RT(1,currentTrial) = TE_behaviour.RTl(1,currentTrial);
        TE_behaviour.FirstLick(1,currentTrial) =  2;
    elseif TE_behaviour.RTr(1,currentTrial) < TE_behaviour.RTl(1,currentTrial)
        TE_behaviour.RT(1,currentTrial) = TE_behaviour.RTr(1,currentTrial);
        TE_behaviour.FirstLick(1,currentTrial) =  1;
    elseif ~isnan(TE_behaviour.RTr(1,currentTrial)) && isnan(TE_behaviour.RTl(1,currentTrial))
        TE_behaviour.RT(1,currentTrial) = TE_behaviour.RTr(1,currentTrial);
        TE_behaviour.FirstLick(1,currentTrial) =  1;
    elseif isnan(TE_behaviour.RTr(1,currentTrial)) && ~isnan(TE_behaviour.RTl(1,currentTrial))
        TE_behaviour.RT(1,currentTrial) = TE_behaviour.RTl(1,currentTrial);
        TE_behaviour.FirstLick(1,currentTrial) =  2;
    end       
end
[TE_behaviour.GoRT, TE_behaviour.NoGoRT]        = deal(TE_behaviour.RT);
TE_behaviour.GoRT(TE_behaviour.TrialType~=1)    = NaN;   % 'reaction time' in the likely rewarded trials
TE_behaviour.NoGoRT(TE_behaviour.TrialType~=2)  = NaN;   % 'reaction time' in the likely punished trials

for currentTrial = 1:ntrials
    switch (TE_behaviour.TrainingStage(1,1))	
      case {1, 2} % TrainingStage 1
        TE_behaviour.Tup{1, currentTrial} = SessionData.RawEvents.Trial{1, currentTrial}.Events.Tup;
        TE_behaviour.StartRewardValveTime(1,currentTrial) = SessionData.RawEvents.Trial{1, currentTrial}.States.Reward(1:1);
        TE_behaviour.EndRewardValveTime(1,currentTrial)   = SessionData.RawEvents.Trial{1, currentTrial}.States.Reward(2);
          
      case 3 % TrainingStage 3
        TE_behaviour.Tup{1, currentTrial} = SessionData.RawEvents.Trial{1, currentTrial}.Events.Tup;
        TE_behaviour.StartRewardSafeValveTime(1,currentTrial) = SessionData.RawEvents.Trial{1, currentTrial}.States.RewardSafe(1:1);
        TE_behaviour.EndRewardSafeValveTime(1,currentTrial)   = SessionData.RawEvents.Trial{1, currentTrial}.States.RewardSafe(2);
        TE_behaviour.StartRewardRiskValveTime(1,currentTrial) = SessionData.RawEvents.Trial{1, currentTrial}.States.RewardRisk(1:1);
        TE_behaviour.EndRewardRiskValveTime(1,currentTrial)   = SessionData.RawEvents.Trial{1, currentTrial}.States.RewardRisk(2);
        TE_behaviour.StartPunishValveTime(1,currentTrial)     = SessionData.RawEvents.Trial{1, currentTrial}.States.Punish(1:1);
        TE_behaviour.EndPunishValveTime(1,currentTrial)       = SessionData.RawEvents.Trial{1, currentTrial}.States.Punish(2);
    end
end

% Sync correction
for currentTrial = 1:ntrials
    if isfield(SessionData.RawEvents.Trial{1,currentTrial}.States,'DeliverStimulus')   % TTL was sent with 25 ms delay by mistake
        TE_behaviour.NeedsSyncCorrection(currentTrial) = 1;
    end
end

% Save
if ifsave == 1
    % For each session, save the individual variables from the extracted TE
    if ~exist([folder filesep sessionID],'dir') == 1
        mkdir([folder filesep sessionID]);
    end
    TEdir = [folder filesep sessionID];
    savename = 'TE_behaviour.mat';
    save([TEdir filesep savename],'-struct','TE_behaviour')
end

disp('TE analysis complete');
end
