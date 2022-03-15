function GamblingTask_reversal
% GAMBLING TASK protocol
%
% -------------------------------------------------------------------------
% Based on CuedOutComeTask.m provided by Heguedus Panna
% and Poking2AFC_3_Test provided by Nicola Solari.
%
% Cecília Pardo-Bellver, 2020
% Laboratory of Network Neurophysiology
% Institute of Experimental Medicine, Hungary.
% 
% Mod. 01/2022 CP-B
% -------------------------------------------------------------------------

% Load Bpod variables
global BpodSystem

% Define parameters -------------------------------------------------------
% Load settings chosen in launch manager into current workspace as a struct called S
S = BpodSystem.ProtocolSettings; %#ok<NASGU> 
S = struct;
if isempty(fieldnames(S))  % If settings file was an empty struct, populate struct with default settings
    
    % Define trinig stage: 1-3 training stages. 1) Lick L & R; 2) Lick to tone; 3) Final.
    list = {'Stage 1', 'Stage 2', 'Stage 3'};
        
    Stg = listdlg('PromptString','Choose Trainig Stage',...
        'SelectionMode','single',...
        'ListString',list);
                       
    switch Stg 
        case 1
            S.TrainingStage = 1;
        case 2
            S.TrainingStage = 2;
        case 3
            S.TrainingStage = 3;
    end
        
    % GUI parameters. It's interactive, you can change it during the task
    S.GUI.TT1_PanicButon  = 0;  % For Trial Type 1 Imposition
    S.GUI.TT2_PanicButon  = 0;  % For Trial Type 2 Imposition
    S.GUI.RewardSmall     = 0.5;  % ul Small reward amount  
    S.GUI.RewardBig       = 8;  % ul Big reward amount 
    S.GUI.SinWavekHz1     = 10; % Cue tone #1 in kHz - tone #1
    S.GUI.SinWavedB1      = 60; % Cue tone #1 dB SPL
    S.GUI.SinWavekHz2     = 04; % Cue tone #2 in kHz - tone #2
    S.GUI.SinWavedB2      = 60; % Cue tone #2 dB SPL
    S.GUI.WhiteWave       = 00; % WhiteNoise tone dB SPL

    % Other Bpod Parameters for the protocol
    S.NoLick        = 1.5;  % in s, the animal needs to learn patience
    S.ITI           = 1;    % ITI duration is set to be exponentially distributed later
    S.SoundDuration = 1;    % in s
    S.RewardValveR  = 1;    % port #1 controls water valve R lickport
    S.RewardValveL  = 4;    % port #2 controls water valve L lickport
    S.PunishValve   = 2;    % port #3 controls air valve
    S.RewardAmountS = S.GUI.RewardSmall; % ul Small reward amount
    S.RewardAmountB = S.GUI.RewardBig;   % ul Big reward amount
    S.PunishTime    = 0.2;  % in s, time with valve open
    S.RewardTimeRS  = GetValveTimes(S.RewardAmountS, 3); % Small reward valve time for right valve
    S.RewardTimeLS  = GetValveTimes(S.RewardAmountS, 1); % Small reward valve time fot left valve
    S.RewardTimeRB  = GetValveTimes(S.RewardAmountB, 3); % Big reward valve time for right valve
    S.RewardTimeLB  = GetValveTimes(S.RewardAmountB, 1); % Big reward valve time for left valve
    S.Type1         = 0.5;          % Probability of trial type 1 (Right)
    S.Type2         = 1 - S.Type1;  % Probability of trial type 2 (Left)
    S.TT1_Panic     = S.GUI.TT1_PanicButon; % For Trial Type 1 Imposition
    S.TT2_Panic     = S.GUI.TT2_PanicButon; % For Trial Type 2 Imposition
    S.SafeProb      = 0.95; % Safe side, probability reward, 5% Omiss ion
    S.RiskProb      = 0.70; % Risky side - 65% big reward , 5% Omission
    S.PunishProb    = 1 - (S.RiskProb + 0.05); % PunishProb = 30% punishment
end
% Initialize parameter GUI plugin
BpodParameterGUI('init', S);

% Define stimuli and send to sound server ---------------------------------
TeensyPort = FindTeensyPort;
T = TeensyAudioPlayer(TeensyPort{1,1}); % Opens Teensy Server at the COM port engaged by the device
FilePath = fullfile(BpodSystem.Path.ProtocolFolder ,  'TeensyCalData.mat');
load(FilePath); % Load Calibration data as a reference
% Generation of Tones 1, 2 and WhiteNiose
AudioPanic  = 0;                    % No WN
SF          = 44100;                % Sound card sampling rate
Tg1         = S.GUI.SinWavedB1;     % Wanted dB for tone 1
Tg2         = S.GUI.SinWavedB2;     % Wanted dB for tone 2
TgWN        = 95;                   % Wanted dB for WN
Fr1         = S.GUI.SinWavekHz1;    % Wanted kHz for tone 1
Fr2         = S.GUI.SinWavekHz2;    % Wanted kHz for tone 2
FrWN        = 22;                   % CHECK! Because the WhiteNoise is saved at the 22nd during the calibration
SPL1        = TeensyCalData.SPL(Fr1);   % Recalls calibrated dB for the frequency of tone 1
SPL2        = TeensyCalData.SPL(Fr2);   % Recalls calibrated dB for the frequency of tone 2
SPLWN       = TeensyCalData.SPL(FrWN);  % Recalls calibrated dB for the WhiteNoise
Ampl1       = TeensyCalData.Amplitude(Fr1);  % Recalls calibrated amplitude for the tone 1 frequency
Ampl2       = TeensyCalData.Amplitude(Fr2);  % Recalls calibrated amplitude for the tone 2 frequency
AmplWN      = TeensyCalData.Amplitude(FrWN); % Recalls calibrated amplitude for the WhiteNoise
NewAmpl1    = AmplAdjst(SPL1,Tg1,Ampl1);    % Calculates new amplitude for tone 1
NewAmpl2    = AmplAdjst(SPL2,Tg2,Ampl2);    % Calculates new amplitude for tone 2
NewAmplWN   = AmplAdjst(SPLWN,TgWN,AmplWN); % Calculates new amplitude for WhiteNoise
sinewave1   = NewAmpl1.*sin(2*pi*Fr1*1000/SF.*(0:SF*S.SoundDuration)); % Creates the sinewave of tone 1
sinewave2   = NewAmpl2.*sin(2*pi*Fr2*1000/SF.*(0:SF*S.SoundDuration)); % Creates the sinewaves of tone 2
sinewaveWN  = NewAmplWN*randn(1,2*44100);                              % Creates the long sinewaves of WhiteNoise
T.load(Fr1, sinewave1);   % Upload to SD card the sinwave
T.load(Fr2, sinewave2);   % Upload to SD card the sinwave
T.load(FrWN, sinewaveWN); % Upload to SD card the sinwave

% Define trials -----------------------------------------------------------
MaxTrials = 300;
rng('shuffle'); % Reset pseudorandom seed

% Initial trials up to 20 are random, but for Stage 1
TrialType  = nan(1,MaxTrials);
TrialType1 = ones(1,(S.Type1*20));      % Generates an array for trial type 1 in the given percentage S.Type1
TrialType2 = repmat(2,1,(S.Type2*20));  % Generates an array for trial type 2 in the given percentage S.Type2
Trials     = [TrialType1, TrialType2];  % Which side should be safe - 1:right, 2:left
TrialType(1:20) = Trials(randperm(length(Trials))); % Random permutation of cued safe side
if S.TrainingStage == 1 % Alternates sides
    TrialType          = ones(1,MaxTrials);
    TrialType(1:2:end) = 2;
elseif S.TrainingStage == 3
    TrialType(1:3)     = 1;
    TrialType(4:6)     = 2;
end
BpodSystem.Data.TrialTypes = []; % The trial type of each trial completed will be added here.

% Defining outcome contingencies, up to 20 are random, but in Stage 1 -----
p = rand([1,length(Trials)]);     % Control outcome contingencies
UsOutcome1 = zeros(size(Trials)); % Safe side outcome
UsOutcome2 = zeros(size(Trials)); % Risky side outcome
if     S.TrainingStage == 3 
    UsOutcome1(p <= S.SafeProb) = 1; % Safe side
    UsOutcome2(p <= S.RiskProb) = 2; % Risk side 
    UsOutcome2(p <= (S.RiskProb + S.PunishProb)& p > S.RiskProb) = 3; % Punishment
    UsOutcome2(1:6) = 2;
elseif S.TrainingStage == 2
    UsOutcome1 = repmat(2,1,length(Trials)); % Small reward always
elseif S.TrainingStage == 1
    UsOutcome1 = repmat(2,1,MaxTrials); % Small reward always
end

% Define delays (between Delaymin and Delaymax in s) ----------------------
rng('shuffle');
if     S.TrainingStage == 2
    Delaymin = 6;
    Delaymax = 9;
elseif S.TrainingStage == 3
    Delaymin = 2;
    Delaymax = 5;
elseif S.TrainingStage == 1
    Delaymin = 0;
    Delaymax = 0;
end
Delay = randi([Delaymin, Delaymax], 1, MaxTrials);
BpodSystem.Data.Delay = Delay;

% Define feedback delay (between Delaymin and Delaymax in s) --------------
if     S.TrainingStage == 2
    FDelaymin = 0.100;
    FDelaymax = 0.300;
elseif S.TrainingStage == 3
    FDelaymin = 0.200;
    FDelaymax = 0.400;
elseif S.TrainingStage == 1
    FDelaymin = 0;
    FDelaymax = 0;
end
FDelay = FDelaymin + (FDelaymax - FDelaymin)*sum(rand(MaxTrials,3),2)/3;
BpodSystem.Data.FDelay = FDelay;

% Initialize plots --------------------------------------------------------
BpodSystem.ProtocolFigures.OutcomePlotFig = figure('Position', [600 600 1100 250],...
    'Name','Outcome plot','numbertitle','off', 'MenuBar', 'none', 'Resize', 'off');
BpodSystem.GUIHandles.OutcomePlot = axes('Position', [0.1300 0.3 0.7750 0.60]);
if S.TrainingStage ==1 || S.TrainingStage ==2
    OutcomePlot_Gambling12(BpodSystem.GUIHandles.OutcomePlot,'init',1-TrialType, UsOutcome1);
else
    OutcomePlot_Gambling(BpodSystem.GUIHandles.OutcomePlot,'init',1-TrialType, [UsOutcome1;UsOutcome2]);
end

% Initialize reward display -----------------------------------------------
TotalRewardDisplay('init');

% Training stage 1 ========================================================
if S.TrainingStage == 1 % Training to lick L & R
          
    for currentTrial = 1:MaxTrials
        
        S = BpodParameterGUI('sync', S); % Sync parameters with BpodParameterGUI plugin
        
        % GUI Update ------------------------------------------------------
        % Controls if parameters for Water reward are changed: if so, it is modified accordingly
        % 1) Check changes in the water rewarded
        newRewardSmall  = S.GUI.RewardSmall; % ul Small reward amount
        if ~isequal(S.RewardAmountS,newRewardSmall)
            S.RewardAmountS = newRewardSmall;
        end
        S.RewardTimeRS  = GetValveTimes(S.RewardAmountS, 3); % Small reward valve time for right valve
        S.RewardTimeLS  = GetValveTimes(S.RewardAmountS, 1); % Small reward valve time fot left valve
        % 2) Sync parameters with BpodParameterGUI plugin
        S = BpodParameterGUI('sync', S);
        
        % Trial conditions ------------------------------------------------
        if     TrialType(currentTrial) == 1
            RewardArg       = S.RewardValveR;
            RewardValveTime = S.RewardTimeRS;
            PortCode        = 'Port1In';
        elseif TrialType(currentTrial) == 2
            RewardArg       = S.RewardValveL;
            RewardValveTime = S.RewardTimeLS;
            PortCode        = 'Port3In';
        end
        
        % Assemble state matrix -------------------------------------------
        sma = NewStateMatrix();
        
        sma = AddState(sma,'Name', 'Start', ...
            'Timer', S.NoLick,...
            'StateChangeConditions', {'Tup', 'WaitforLick','Port1In','Start','Port3In','Start'},...
            'OutputActions', {});
        
        sma = AddState(sma,'Name', 'Reward', ...
            'Timer',RewardValveTime,...
            'StateChangeConditions', {'Tup', 'PostUS'},...
            'OutputActions', {'ValveState', RewardArg,'BNC2', 1});  
       
        sma = AddState(sma,'Name','WaitforLick',...
            'Timer',0,...
            'StateChangeConditions',{PortCode,'FeedBackDelay'},...
            'OutputActions',{});
        
        sma = AddState(sma,'Name','FeedBackDelay',...
            'Timer',FDelay(currentTrial),...
            'StateChangeConditions',{'Tup','Reward'},...
            'OutputActions',{});
        
        sma = AddState(sma,'Name','PostUS',...
            'Timer',1,...
            'StateChangeConditions',{'Port1In','ResetDrinkingTimer', 'Port3In','ResetDrinkingTimer', 'Tup','exit'},...
            'OutputActions',{'BNC2', 0}); % 
        
        sma = AddState(sma,'Name','ResetDrinkingTimer',...
            'Timer',0,...
            'StateChangeConditions',{'Tup','PostUS'},...
            'OutputActions',{}); % keep the animal in PostUS until licking stops for 1 s
        SendStateMatrix(sma);
        RawEvents = RunStateMatrix;     
        
        % Update trial data -----------------------------------------------
        if ~isempty(fieldnames(RawEvents)) % If trial data was returned
            % Computes trial events from raw data
            BpodSystem.Data = AddTrialEvents(BpodSystem.Data,RawEvents); 
            % Adds the settings used for the current trial to the Data struct
            BpodSystem.Data.TrialSettings(currentTrial) = S; 
            % Adds the trial type of the current trial to data
            BpodSystem.Data.TrialTypes(currentTrial) = TrialType(currentTrial); 
            % Update outcome plot
            Outcomes = UpdateOutcomePlot(TrialType, BpodSystem.Data, UsOutcome1);
            % Adds the outcome of the current trial to data
            BpodSystem.Data.TrialOutcome(currentTrial) = Outcomes(currentTrial);
            % Saves the field BpodSystem.Data to the current data file
            SaveBpodSessionData;
        end
        
        HandlePauseCondition;
        if BpodSystem.Status.BeingUsed == 0
            return    
        end
        
        % Display water consumed ------------------------------------------
        SafeChoice  = nan(1,BpodSystem.Data.nTrials);
        for ii = 1:BpodSystem.Data.nTrials
            if     Outcomes(ii) == 1
                SafeChoice(ii) = 1;
            elseif Outcomes(ii) == 2
                SafeChoice(ii) =  1;
            else
            end
        end
        if SafeChoice(currentTrial) == 1
            WaterReceived = S.RewardAmountS;
            TotalRewardDisplay('add', WaterReceived);
        elseif isnan(SafeChoice(currentTrial)) == 1
            WaterReceived = 0;
            TotalRewardDisplay('add', WaterReceived);
        end
    end
end

% Training stage 2 ========================================================
if S.TrainingStage == 2 %introducing a cue and delay
    
    for currentTrial = 1:MaxTrials
        
        S = BpodParameterGUI('sync', S); % Synchronize the GUI
        
        % GUI Update ------------------------------------------------------
        % Controls if parameters for Tone/ Water reward are changed: if so, it is modified accordingly
        % 1) Check changes in the audio
        newTg1 = S.GUI.SinWavedB1;
        newTg2 = S.GUI.SinWavedB2;
        newFr1 = S.GUI.SinWavekHz1;
        newFr2 = S.GUI.SinWavekHz2;
        if     ~isequal(Tg1,newTg1) || ~isequal(Fr1,newFr1) 
            SPL1      = TeensyCalData.SPL(S.GUI.SinWavekHz1);
            Ampl1     = TeensyCalData.Amplitude(S.GUI.SinWavekHz1);
            NewAmpl1  = AmplAdjst(SPL1,newTg1,Ampl1);
            sinewave1 = NewAmpl1.*sin(2*pi*newFr1*1000/SF.*(0:SF*S.SoundDuration));
            TeensySoundServer('loadwaveform', S.GUI.SinWavekHz1, sinewave1);
            Tg1       = newTg1;
            Fr1       = newFr1;
        elseif ~isequal(Tg2,newTg2) || ~isequal(Fr2,newFr2) 
            SPL2      = TeensyCalData.SPL(S.GUI.SinWavekHz2);
            Ampl2     = TeensyCalData.Amplitude(S.GUI.SinWavekHz2);
            NewAmpl2  = AmplAdjst(SPL2,newTg2,Ampl2);
            sinewave2 = NewAmpl2.*sin(2*pi*newFr2*1000/SF.*(0:SF*S.SoundDuration));
            TeensySoundServer('loadwaveform', S.GUI.SinWavekHz2, sinewave2);
            Tg2       = newTg2;
            Fr2       = newFr2;
        end
        % 2) Check changes in the water rewarded
        newRewardSmall  = S.GUI.RewardSmall; % ul Small reward amount
        newRewardBig    = S.GUI.RewardBig; % ul Big reward amount
        if ~isequal(S.RewardAmountS,newRewardSmall) || ~isequal(S.RewardAmountB,newRewardBig)
            S.RewardAmountS = newRewardSmall;
            S.RewardAmountB = newRewardBig;
        end
        S.RewardTimeRS  = GetValveTimes(S.RewardAmountS, 3); % Small reward valve time for right valve
        S.RewardTimeLS  = GetValveTimes(S.RewardAmountS, 1); % Small reward valve time fot left valve
        % 3) Check TrialType imposition
        TT1_PanicButon = S.GUI.TT1_PanicButon; % For Trial Type 1 Imposition
        TT2_PanicButon = S.GUI.TT2_PanicButon; % For Trial Type 2 Imposition
        if TT1_PanicButon == 1 || TT2_PanicButon == 1
            S.TT1_Panic = TT1_PanicButon;
            S.TT2_Panic = TT2_PanicButon;
        end
        % 4) Check WN panic button
        WNPanic = S.GUI.WhiteWave;
        if WNPanic == 1
            AudioPanic = 1;
        else 
        end
        % 5) Sync parameters with BpodParameterGUI plugin
        S = BpodParameterGUI('sync', S); 
        
        % Trial Type decision ---------------------------------------------
        if currentTrial <= 20 % first 20 trials are purely random
            Correction     = 0;
            TT1_PanicButon = 0;
            TT2_PanicButon = 0;
        else
            Correction = BiasCounter(TrialType, BpodSystem.Data);
        end
        
        if currentTrial > 20
            LastTrials = sum(TrialType(currentTrial-5:currentTrial-1));
            if     TT1_PanicButon == 1 || Correction == 1
                TrialType(currentTrial) = 1;
            elseif TT2_PanicButon == 1 || Correction == 2
                TrialType(currentTrial) = 2;
            elseif LastTrials == 5
                TrialType(currentTrial) = 2;
            elseif LastTrials == 10
                TrialType(currentTrial) = 1;
            else 
                TrialType(currentTrial) = ceil(rand(1,1)*2);
            end
        else
        end
            
        % Current trial conditions ----------------------------------------
        if     TrialType(currentTrial) == 2
            Audio           = Fr1;
            PortCondition   = 'Port3In';
            RewardValveTime = S.RewardTimeLS;
            RewardValveCode = S.RewardValveL;
        elseif TrialType(currentTrial) == 1
            Audio           = Fr2;
            PortCondition   = 'Port1In';
            RewardValveTime = S.RewardTimeRS;
            RewardValveCode = S.RewardValveR;
        end
        if AudioPanic == 1
            Audio = FrWN;
        end
        
        % Inter-trial interval --------------------------------------------
        S.ITI = 10;
        while S.ITI > 3   % ITI dustribution: 1 + exponential, truncated at 4 s
            S.ITI = exprnd(1)+1;
        end
        
        % Assemble state matrix -------------------------------------------
        sma = NewStateMatrix();
        
        sma = AddState(sma,'Name', 'NoLick', ...
            'Timer', S.NoLick,...
            'StateChangeConditions', {'Tup', 'ITI','Port1In','RestartNoLick', 'Port3In','RestartNoLick'},...
            'OutputActions', {});
        
        sma = AddState(sma,'Name', 'RestartNoLick', ...
            'Timer', 0,...
            'StateChangeConditions', {'Tup', 'NoLick',},...
            'OutputActions', {});
        
        sma = AddState(sma, 'Name', 'ITI', ...
            'Timer',S.ITI,...
            'StateChangeConditions', {'Tup', 'StartStimulus','Port1In','RestartNoLick', 'Port3In','RestartNoLick'},...
            'OutputActions', {});  
        
        sma = AddState(sma, 'Name', 'StartStimulus', ...
            'Timer', S.SoundDuration,...
            'StateChangeConditions', {PortCondition,'FeedBackDelay','Tup','Delay'},...
            'OutputActions', {'TeensyAudio1', Audio, 'BNC2', 1});   
        
        sma = AddState(sma,'Name','FeedBackDelay',...
            'Timer',FDelay(currentTrial),...
            'StateChangeConditions',{'Tup','Reward'},...
            'OutputActions',{});
        
        sma = AddState(sma, 'Name','Delay', ...
            'Timer', Delay(currentTrial),...
            'StateChangeConditions', {PortCondition,'FeedBackDelay','Tup','PostUS'},...
            'OutputActions', {});
        
        sma = AddState(sma,'Name', 'Reward', ...
            'Timer',RewardValveTime,...
            'StateChangeConditions', {'Tup', 'PostUS'},...
            'OutputActions', {'ValveState', RewardValveCode, 'BNC2', 0});  
        
        sma = AddState(sma,'Name','PostUS',...
            'Timer',1,...
            'StateChangeConditions',{PortCondition,'ResetDrinkingTimer', 'Tup','exit'},...
            'OutputActions',{'BNC2', 0});   
        
        sma = AddState(sma,'Name','ResetDrinkingTimer',...
            'Timer',0,...
            'StateChangeConditions',{'Tup','PostUS'},...
            'OutputActions',{});   % keep the animal in PostUS until licking stops for 1 s
        SendStateMatrix(sma);
        RawEvents = RunStateMatrix;
        
        % Update trial data -----------------------------------------------
        if ~isempty(fieldnames(RawEvents)) % If trial data was returned
            % Computes trial events from raw data
            BpodSystem.Data = AddTrialEvents(BpodSystem.Data,RawEvents); 
            % Adds the settings used for the current trial to the Data struct
            BpodSystem.Data.TrialSettings(currentTrial) = S; 
            % Adds the trial type of the current trial to data
            BpodSystem.Data.TrialTypes(currentTrial) = TrialType(currentTrial); 
            % Update outcome plot
            Outcomes = UpdateOutcomePlot(TrialType, BpodSystem.Data, UsOutcome1); 
            % Adds the outcome of the current trial to data
            BpodSystem.Data.TrialOutcome(currentTrial) = Outcomes(currentTrial);
            % Saves the field BpodSystem.Data to the current data file 
            SaveBpodSessionData;    
        end
        
        HandlePauseCondition;
        if BpodSystem.Status.BeingUsed == 0
            return
        end
        
        % Display water consumed ------------------------------------------
        SafeChoice  = nan(1,BpodSystem.Data.nTrials);
        for ii = 1:BpodSystem.Data.nTrials
            if     Outcomes(ii) == 1
                SafeChoice(ii) = 1;
            elseif Outcomes(ii) == 2
                SafeChoice(ii) =  1;
            else
            end
        end
        if SafeChoice(currentTrial) == 1
            WaterReceived = S.RewardAmountS;
            TotalRewardDisplay('add', WaterReceived);
        elseif isnan(SafeChoice(currentTrial)) == 1
            WaterReceived = 0;
            TotalRewardDisplay('add', WaterReceived);
        end
    end    
end

% Training stage 3 ========================================================
if S.TrainingStage == 3 %introducing a cue and delay
    
    for currentTrial = 1:MaxTrials
        
        S = BpodParameterGUI('sync', S); % Synchronize the GUI
        
        % GUI Update ------------------------------------------------------
        % Controls if parameters for Tone/ WaterReward are changed: if so, it is modified
        % 1) Check cahnges in the audio
        newTg1 = S.GUI.SinWavedB1;
        newTg2 = S.GUI.SinWavedB2;
        newFr1 = S.GUI.SinWavekHz1;
        newFr2 = S.GUI.SinWavekHz2;
        if     ~isequal(Tg1,newTg1) || ~isequal(Fr1,newFr1) 
            SPL1      = TeensyCalData.SPL(S.GUI.SinWavekHz1);
            Ampl1     = TeensyCalData.Amplitude(S.GUI.SinWavekHz1);
            NewAmpl1  = AmplAdjst(SPL1,newTg1,Ampl1);
            sinewave1 = NewAmpl1.*sin(2*pi*newFr1*1000/SF.*(0:SF*S.SoundDuration));
            TeensySoundServer('loadwaveform', S.GUI.SinWavekHz1, sinewave1);
            Tg1       = newTg1;
            Fr1       = newFr1;
        elseif ~isequal(Tg2,newTg2) || ~isequal(Fr2,newFr2) 
            SPL2      = TeensyCalData.SPL(S.GUI.SinWavekHz2);
            Ampl2     = TeensyCalData.Amplitude(S.GUI.SinWavekHz2);
            NewAmpl2  = AmplAdjst(SPL2,newTg2,Ampl2);
            sinewave2 = NewAmpl2.*sin(2*pi*newFr2*1000/SF.*(0:SF*S.SoundDuration));
            TeensySoundServer('loadwaveform', S.GUI.SinWavekHz2, sinewave2);
            Tg2       = newTg2;
            Fr2       = newFr2;
        end
        % 2) Check changes in the water rewarded
        newRewardSmall  = S.GUI.RewardSmall; % ul Small reward amount
        newRewardBig    = S.GUI.RewardBig; % ul Big reward amount
        if ~isequal(S.RewardAmountS,newRewardSmall) || ~isequal(S.RewardAmountB,newRewardBig)
            S.RewardAmountS = newRewardSmall;
            S.RewardAmountB = newRewardBig;
        end
        S.RewardTimeRS  = GetValveTimes(S.RewardAmountS, 3); % Small reward valve time for right valve
        S.RewardTimeLS  = GetValveTimes(S.RewardAmountS, 1); % Small reward valve time fot left valve
        % 3) Check TrialType imposition
        TT1_PanicButon = S.GUI.TT1_PanicButon; % For Trial Type 1 Imposition
        TT2_PanicButon = S.GUI.TT2_PanicButon; % For Trial Type 2 Imposition
        if TT1_PanicButon == 1 || TT2_PanicButon == 1
            S.TT1_Panic = TT1_PanicButon;
            S.TT2_Panic = TT2_PanicButon;
        end
        % 4) Check WN panic button
        WNPanic = S.GUI.WhiteWave;
        if WNPanic == 1
            AudioPanic = 1;
        else 
        end
        % 5) Sync parameters with BpodParameterGUI plugin
        S = BpodParameterGUI('sync', S); 
        
        % Trial Type decision ---------------------------------------------
        if currentTrial <= 20 % first 20 trials are purely random
            Correction     = 0;
            TT1_PanicButon = 0;
            TT2_PanicButon = 0;
        else
            Correction = BiasCounter(TrialType, BpodSystem.Data);
        end
        
        if currentTrial > 20
            LastTrials = sum(TrialType(currentTrial-5:currentTrial-1));
            if     TT1_PanicButon == 1 
                TrialType(currentTrial) = 1;
            elseif TT2_PanicButon == 1 
                TrialType(currentTrial) = 2;
            elseif Correction == 1
                TrialType(currentTrial) = 1;
            elseif Correction == 2
                TrialType(currentTrial) = 2;
            elseif LastTrials == 5
                TrialType(currentTrial) = 2;
            elseif LastTrials == 10
                TrialType(currentTrial) = 1;
            else
                TrialType(currentTrial) = ceil(rand(1,1)*2);
            end
        else 
        end
        
        % UsOutcome decision ----------------------------------------------
        clearvars A B p
        if currentTrial/5 == round(currentTrial/5)
            rng('shuffle');
        end
        A = 0;
        B = 0;
        if currentTrial > 20
            p = rand(1,1);
            if UsOutcome1(currentTrial-1) == 0
                A = 1;
            elseif p >= (1- S.SafeProb)
                A = 1;
            end
            if UsOutcome2(currentTrial-1) == 0
                B = 2;             
            elseif p <= (S.RiskProb + S.PunishProb) && p > S.RiskProb
                B = 3;
            elseif p <= S.RiskProb
                B = 2;
            end
            UsOutcome1(currentTrial) = A; % Safe side outcome
            UsOutcome2(currentTrial) = B; % Risk side outcome
        else
        end
        
        % Current trial conditions ----------------------------------------
        if     TrialType(currentTrial) == 2 %choosing port
            Audio               = Fr1;
            SafePort            = 'Port1In';
            RiskyPort           = 'Port3In';
            RewardValveTimeSafe = S.RewardTimeRS; %safe right port small reward size
            RewardValveTimeRisk = S.RewardTimeLB; %risk left port big reward size
            RewardValveCodeSafe = S.RewardValveR;
            RewardValveCodeRisk = S.RewardValveL;
        elseif TrialType(currentTrial) == 1
            Audio               = Fr2;
            SafePort            = 'Port3In';
            RiskyPort           = 'Port1In';
            RewardValveTimeRisk = S.RewardTimeRB; %risk right port big reward size
            RewardValveTimeSafe = S.RewardTimeLS; %safe left port small reward size
            RewardValveCodeSafe = S.RewardValveL;
            RewardValveCodeRisk = S.RewardValveR;
        end
        if AudioPanic == 1
            Audio = FrWN;
        end
        
        % Adjust Safe/ Risk according to Outcomes -------------------------
        if     UsOutcome1(currentTrial) == 1
            SafePortArgument = 'RewardSafe'; 
        elseif UsOutcome1(currentTrial) == 0
            SafePortArgument = 'PostUS';
        end
        
        if     UsOutcome2(currentTrial) == 2
            RiskyPortArgument = 'RewardRisk';
        elseif UsOutcome2(currentTrial) == 3
            RiskyPortArgument = 'Punish';
        elseif UsOutcome2(currentTrial) == 0
            RiskyPortArgument = 'PostUS';
        end
        
        % Inter-trial interval --------------------------------------------
        S.ITI = 10;
        while S.ITI > 4   % ITI dustribution: 1 + exponential, truncated at 4 s
            S.ITI = exprnd(1)+1;
        end
        
        % Assemble state matrix -------------------------------------------
        sma = NewStateMatrix();
        
        sma = AddState(sma,'Name', 'NoLick', ...
            'Timer', S.NoLick,...
            'StateChangeConditions', {'Tup', 'ITI','Port1In','RestartNoLick','Port3In','RestartNoLick'},...
            'OutputActions', {}); % Animal in no lick state until S.NoLick pause in licking & no licks on Port1
        
        sma = AddState(sma,'Name', 'RestartNoLick', ...
            'Timer', 0,...
            'StateChangeConditions', {'Tup', 'NoLick'},...
            'OutputActions', {}); % Resets 'NoLick' state to 0 if animal licks
        
        sma = AddState(sma, 'Name', 'ITI', ...
            'Timer',S.ITI,...
            'StateChangeConditions', {'Tup', 'StartStimulus', 'Port1In','RestartNoLick', 'Port3In','RestartNoLick'},...
            'OutputActions', {});  % 1 + exponential foreperiod < 4s
        
        sma = AddState(sma, 'Name', 'StartStimulus', ...
            'Timer', S.SoundDuration,...
            'StateChangeConditions', {SafePort,'SafeFeedBackDelay',RiskyPort,'RiskFeedBackDelay','Tup','Delay'},...
            'OutputActions', {'TeensyAudio1', Audio, 'BNC2', 1});
        
        sma = AddState(sma,'Name','SafeFeedBackDelay',...
            'Timer',FDelay(currentTrial),...
            'StateChangeConditions',{'Tup',SafePortArgument},...
            'OutputActions',{});
        
        sma = AddState(sma,'Name','RiskFeedBackDelay',...
            'Timer',FDelay(currentTrial),...
            'StateChangeConditions',{'Tup',RiskyPortArgument},...
            'OutputActions',{});
        
        sma = AddState(sma, 'Name','Delay', ...
            'Timer', Delay(currentTrial),...
            'StateChangeConditions', {SafePort,SafePortArgument, RiskyPort,RiskyPortArgument,'Tup','PostUS'},...
            'OutputActions', {});
        
        sma = AddState(sma,'Name', 'RewardSafe', ...
            'Timer',RewardValveTimeSafe,...
            'StateChangeConditions', {'Tup', 'PostUS'},...
            'OutputActions', {'ValveState', RewardValveCodeSafe,'BNC2', 0});   % deliver water
        
        sma = AddState(sma,'Name', 'RewardRisk', ...
            'Timer',RewardValveTimeRisk,...
            'StateChangeConditions', {'Tup', 'PostUS'},...
            'OutputActions', {'ValveState', RewardValveCodeRisk,'BNC2', 0});   % deliver water
        
        sma = AddState(sma, 'Name', 'Punish', ...
            'Timer',S.PunishTime, ...
            'StateChangeConditions', {'Tup', 'PostUS'}, ...
            'OutputActions', {'ValveState', S.PunishValve,'BNC2', 0}); %deiliver air puff
        
        sma = AddState(sma,'Name','PostUS',...
            'Timer',1,...
            'StateChangeConditions',{SafePort,'ResetDrinkingTimer', RiskyPort,'ResetDrinkingTimer','Tup','exit'},...
            'OutputActions',{'BNC2', 0});   % drinking
        
        sma = AddState(sma,'Name','ResetDrinkingTimer',...
            'Timer',0,...
            'StateChangeConditions',{'Tup','PostUS'},...
            'OutputActions',{});   % keep the animal in PostUS until licking stops for 1 s
        SendStateMatrix(sma);
        RawEvents = RunStateMatrix;
        
        % Update trial data -----------------------------------------------
        if ~isempty(fieldnames(RawEvents)) % If trial data was returned
            % Computes trial events from raw data
            BpodSystem.Data = AddTrialEvents(BpodSystem.Data,RawEvents); 
            % Adds the settings used for the current trial to the Data struct
            BpodSystem.Data.TrialSettings(currentTrial) = S; 
            % Adds the trial type of the current trial to data
            BpodSystem.Data.TrialTypes(currentTrial) = TrialType(currentTrial); 
            % Update outcome plot
            Outcomes = UpdateOutcomePlot2(TrialType, BpodSystem.Data, UsOutcome1, UsOutcome2);
            % Adds the outcome of the current trial to data
            BpodSystem.Data.TrialOutcome(currentTrial) = Outcomes(currentTrial);
            % Saves the field BpodSystem.Data to the current data file
            SaveBpodSessionData; 
        end
        
        HandlePauseCondition;
        if BpodSystem.Status.BeingUsed == 0
            return    
        end
        
        % Display water consumed ------------------------------------------
        SafeChoice  = nan(1,BpodSystem.Data.nTrials);
        RiskyChoice = nan(1,BpodSystem.Data.nTrials);
        for ii = 1:BpodSystem.Data.nTrials
            if     Outcomes(ii) == 1
                SafeChoice(ii)  = 1;
            elseif Outcomes(ii) == 2
                RiskyChoice(ii) =  1;
            else
            end
        end
        if     SafeChoice(currentTrial) == 1
            WaterReceived = S.RewardAmountS;
            TotalRewardDisplay('add', WaterReceived);
        elseif RiskyChoice(currentTrial) == 1
            WaterReceived = S.RewardAmountB;
            TotalRewardDisplay('add', WaterReceived);
        elseif isnan(SafeChoice(currentTrial)) == 1 || isnan(RiskChoice(currentTrial)) == 1
            WaterReceived = 0;
            TotalRewardDisplay('add', WaterReceived);
        end
    end
end
end

% -------------------------------------------------------------------------
function Outcomes = UpdateOutcomePlot(TrialType, Data, UsOutcome)

% Load Bpod variables
global BpodSystem
Outcomes = zeros(1,Data.nTrials);

% Outcome
for x = 1:Data.nTrials
    if ~isnan(Data.RawEvents.Trial{x}.States.Reward(1))
        Outcomes(x) = 2;   % lick, reward
    else
        Outcomes(x) = 0;   % omission
    end
end
OutcomePlot_Gambling12(BpodSystem.GUIHandles.OutcomePlot,'update',...
    Data.nTrials+1,1-TrialType,Outcomes, UsOutcome)
end

% -------------------------------------------------------------------------
function Outcomes = UpdateOutcomePlot2(TrialType, Data, UsOutcome1, UsOutcome2)

% Load Bpod variablesclc
global BpodSystem
Outcomes = zeros(1,Data.nTrials);
Choice = zeros(1,Data.nTrials);
% Outcome
for x = 1:Data.nTrials
    if     UsOutcome1(x) == 1  && ~isnan(Data.RawEvents.Trial{x}.States.RewardSafe(1))
        Outcomes(x) = 1;   % lick, safe reward
    elseif UsOutcome2(x) == 2  && ~isnan(Data.RawEvents.Trial{x}.States.RewardRisk(1))
        Outcomes(x) = 2;   % lick, risky reward 
    elseif UsOutcome2(x) == 3  && ~isnan(Data.RawEvents.Trial{x}.States.Punish(1)) 
        Outcomes(x) = 3;   % lick, punish (risky port)
    elseif UsOutcome2(x) == 3  && isnan(Data.RawEvents.Trial{x}.States.Punish(1)) 
        Outcomes(x) = 5;   % no lick, punish 
    elseif UsOutcome1(x) == 0  || UsOutcome2(x) == 0
        StartStim = Data.RawEvents.Trial{x}.States.StartStimulus(1);
        EndDelay  = Data.RawEvents.Trial{x}.States.Delay(2);
        Port1 = 0;
        Port3 = 0;
        
        try % Port 1 = R
            Port1A = Data.RawEvents.Trial{x}.Events.Port1In > StartStim;
            Port1B = Data.RawEvents.Trial{x}.Events.Port1In < EndDelay;
            Port1  = sum(Port1A + Port1B);
        catch
        end
        
        try % Port 3 = L
            Port3A = Data.RawEvents.Trial{x}.Events.Port3In > StartStim;
            Port3B = Data.RawEvents.Trial{x}.Events.Port3In < EndDelay;
            Port3  = sum(Port3A + Port3B);
        catch
        end
        
        if     Port3 > 0  || Port1 > 0
            Outcomes(x) = 4;  % lick, omission
        elseif Port3 == 0  && Port1 == 0
            Outcomes(x) = 6;  % no lick, omission
        end
    elseif UsOutcome1(x) == 1  && isnan(Data.RawEvents.Trial{x}.States.RewardSafe(1))
        Outcomes(x) = 0;   % no lick, no reward
    elseif UsOutcome2(x) == 2  && isnan(Data.RawEvents.Trial{x}.States.RewardRisk(1))
        Outcomes(x) = 0;   % no lick, no reward
    end
    
    if     ~isnan(Data.RawEvents.Trial{x}.States.RewardSafe(1))
        Choice(x) = 1;
    elseif ~isnan(Data.RawEvents.Trial{x}.States.RewardRisk(1))
        Choice(x) = 2;
    elseif ~isnan(Data.RawEvents.Trial{x}.States.Punish(1))
        Choice(x) = 2;
    else
        Choice(x) = 1.5;
    end
    
end
OutcomePlot_Gambling(BpodSystem.GUIHandles.OutcomePlot,'update',...
    Data.nTrials+1,1-TrialType, 1-Choice ,Outcomes, [UsOutcome1;UsOutcome2])
end

% -------------------------------------------------------------------------
function [Correction] = BiasCounter(TrialType, Data) % mild correction for bias (...maybe.)

global BpodSystem %#ok<NUSED>
window = 19;
PortChoice = repmat(9,1,Data.nTrials);

for jj = 1:Data.nTrials
    StartStim = Data.RawEvents.Trial{jj}.States.StartStimulus(1);
    EndDelay  = Data.RawEvents.Trial{jj}.States.Delay(2);
        
    try % Right licks - TT1
        Port1A = Data.RawEvents.Trial{jj}.Events.Port1In > StartStim;
        Port1B = Data.RawEvents.Trial{jj}.Events.Port1In < EndDelay;
        if ~isnan(sum(Port1A + Port1B))
            PortChoice(jj) = 1;
        end
    catch
    end

    try % Left lick - TT2
        Port3A = Data.RawEvents.Trial{jj}.Events.Port3In > StartStim;
        Port3B = Data.RawEvents.Trial{jj}.Events.Port3In < EndDelay;
        if ~isnan(sum(Port3A + Port3B))
            PortChoice(jj) = 3;
        end
    catch
    end
end
LastOutcomes = PortChoice(end-(window-1):end);
t1 = tabulate(LastOutcomes);

LastTrials = TrialType(Data.nTrials-(window-1):Data.nTrials);
t2 = tabulate(LastTrials);

if     t2(1,3) > 90 % >90% TT1
    C = 2;
elseif t2(2,3) > 90 % >90% TT2
    C = 1;
elseif t1(1,3) > 60 % >60% Licks R
    C = 2;
elseif t1(2,3) > 60 % >60% Licks L
    C = 1;
else
    C = 0;
end

Correction = C;

end

% -------------------------------------------------------------------------
function [CalAmpl] = AmplAdjst(SPL,Tg,Ampl) % Calculate the new proper sinewave amplitude

y = SPL - Tg;
b =  20 * log10(Ampl) - y;
c = b / 20;
CalAmpl = 10 .^ c;
end
