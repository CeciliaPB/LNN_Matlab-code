function GamblingTask
% GAMBLING TASK protocol
%
% -------------------------------------------------------------------------
% Based on CuedOutComeTask.m provided by Heguedus Panna
% and Poking2AFC_3_Test provided by Nicola Solari.
%
% Cecília Pardo-Bellver, 2020
% Laboratory of Network Neurophysiology
% Institute of Experimantal Medicine, Hungary.
% -------------------------------------------------------------------------

% Load Bpod variables
global BpodSystem

% Define parameters -------------------------------------------------------
S = BpodSystem.ProtocolSettings; % Load settings chosen in launch manager into current workspace as a struct called S
S = struct;
if isempty(fieldnames(S))  % If settings file was an empty struct, populate struct with default settings
    
    S.TrainingStage = 3;    % 1-3 training stages. 1) Lick L & R; 2) Lick to tone; 3) Final.
    
    % GUI parameters. It's interactive, you can change it during the task
    S.GUI.TT1_PanicButon  = 0;  % For Trial Type 1 Imposition
    S.GUI.TT2_PanicButon  = 0;  % For Trial Type 2 Imposition
    S.GUI.RewardSmall     = 2;  % ul Small reward amount  
    S.GUI.RewardBig       = 8;  % ul Big reward amount 
    S.GUI.SinWavekHz1     = 10; % Cue tone #1 in kHz - tone #1
    S.GUI.SinWavedB1      = 60; % Cue tone #1 dB SPL
    S.GUI.SinWavekHz2     = 04; % Cue tone #2 in kHz - tone #2
    S.GUI.SinWavedB2      = 60; % Cue tone #2 dB SPL
    S.GUI.WhiteWavedB     = 95; % WhiteNoise tone dB SPL

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
    S.NoPunish      = false;% true = no airpuff
    S.Type1         = 0.5;  % Probability of trial type 1 (Right)
    S.Type2         = 0.5;  % Probability of trial type 2 (Left)
    S.TT1_Panic     = S.GUI.TT1_PanicButon; % For Trial Type 1 Imposition
    S.TT2_Panic     = S.GUI.TT2_PanicButon; % For Trial Type 2 Imposition
    S.SafeProb      = 0.95; % Safe side, probability reward
    S.RiskProb      = 0.60; % Risky side - 60% big reward
    S.PunishProb    = 0.70; % PunishProb = 100 - 70 = 30% punishment
end
% Initialize parameter GUI plugin
BpodParameterGUI('init', S);

% Define stimuli and send to sound server ---------------------------------
T = TeensyAudioPlayer('COM1'); % Opens Teensy Server at the COM port engaged by the device
FilePath = fullfile(BpodSystem.Path.ProtocolFolder ,  'TeensyCalData.mat');
load(FilePath); % Load Calibration data as a reference
% Generation of Tones 1, 2 and WhiteNiose
SF          = 44100;                % Sound card sampling rate
Tg1         = S.GUI.SinWavedB1;     % Wanted dB for tone 1
Tg2         = S.GUI.SinWavedB2;     % Wanted dB for tone 2
TgWN        = S.GUI.WhiteWavedB;    % Wanted dB for WN
Fr1         = S.GUI.SinWavekHz1;    % Wanted kHz for tone 1
Fr2         = S.GUI.SinWavekHz2;    % Wanted kHz for tone 2
FrWN        = 22;                   % CHECK! Because the WhiteNoise is saved at the 22nd during the calibration (for me, having done it for 21 tracks)
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
% T.load(Fr1, sinewave1);   % Upload to SD card the sinwave
% T.load(Fr2, sinewave2);   % Upload to SD card the sinwave
% T.load(FrWN, sinewaveWN); % Upload to SD card the sinwave

% Define trials -----------------------------------------------------------
MaxTrials = 1000;
rng('shuffle'); % Reset pseudorandom seed

% Initial trials up to 20 are random, but for Stage 1
TrialType  = nan(1,MaxTrials);
TrialType1 = ones(1,(S.Type1*20));      % Generates an array for trial type 1 in the given percentage S.Type1
TrialType2 = repmat(2,1,(S.Type2*20));  % Generates an array for trial type 2 in the given percentage S.Type2
Trials     = [TrialType1, TrialType2];  % Which side should be safe - 1:right, 2:left
TrialType(1:20) = Trials(randperm(length(Trials))); % Random permutation of cued safe side
if S.TrainingStage == 1 
    TrialType          = ones(1,MaxTrials);
    TrialType(1:2:end) = 2;
end
BpodSystem.Data.TrialTypes = []; % The trial type of each trial completed will be added here.

% Defining outcome contingencies, up to 20 are random, but Stage 1 --------
p = rand(1,length(Trials));       % Control outcome contingencies
UsOutcome1 = zeros(size(Trials)); % Safe side outcome
UsOutcome2 = zeros(size(Trials)); % Risky side outcome
if     S.TrainingStage == 3 
    UsOutcome1(p <= S.SafeProb) = 1;            % Safe side
    UsOutcome2(p <= S.RiskProb) = 2;            % Risk side - 60% big reward
    UsOutcome2(p > S.PunishProb & p <= 1) = 3;  % 30% punishment
elseif S.TrainingStage == 2
    UsOutcome1 = repmat(2,1,length(Trials)); % Small reward always
elseif S.TrainingStage == 1
    UsOutcome1 = repmat(2,1,MaxTrials); % Small reward always
end

% Define delays (between Delaymin and Delaymax in s) ----------------------
    if     S.TrainingStage == 2
        Delaymin = 6; 
        Delaymax = 9; 
    elseif S.TrainingStage == 3
        Delaymin = 2; 
        Delaymax = 5; 
    elseif S.TrainingStage == 1
        Delaymin = 1; 
        Delaymax = 2; 
    end
Delay = randi([Delaymin, Delaymax], 1, MaxTrials );
BpodSystem.Data.Delay = Delay;

% Initialize plots --------------------------------------------------------
BpodSystem.ProtocolFigures.OutcomePlotFig = figure('Position', [600 600 1100 250],'Name','Outcome plot','numbertitle','off', 'MenuBar', 'none', 'Resize', 'off');
BpodSystem.GUIHandles.OutcomePlot = axes('Position', [0.1300 0.3 0.7750 0.60]);
if S.TrainingStage ==1 || S.TrainingStage ==2
    OutcomePlot_Gambling12(BpodSystem.GUIHandles.OutcomePlot,'init',1-TrialType, UsOutcome1);
else
    OutcomePlot_Gambling(BpodSystem.GUIHandles.OutcomePlot,'init',1-TrialType, [UsOutcome1;UsOutcome2]);
end

% Training stage 1 --------------------------------------------------------
if S.TrainingStage == 1 % Training to lick L & R
    
    S = BpodParameterGUI('sync', S); % Sync parameters with BpodParameterGUI plugin
    
    for currentTrial = 1:MaxTrials
        
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
            'OutputActions', {'ValveState', RewardArg, 'BNC1', 1});   % deliver water
       
        sma = AddState(sma,'Name','WaitforLick',...
            'Timer',0,...
            'StateChangeConditions',{PortCode,'Reward'},...
            'OutputActions',{});
        
        sma = AddState(sma,'Name','PostUS',...
            'Timer',1,...
            'StateChangeConditions',{'Port1In','ResetDrinkingTimer', 'Port3In','ResetDrinkingTimer', 'Tup','exit'},...
            'OutputActions',{'BNC1', 0,}); % drinking
        
        sma = AddState(sma,'Name','ResetDrinkingTimer',...
            'Timer',0,...
            'StateChangeConditions',{'Tup','PostUS'},...
            'OutputActions',{}); % keep the animal in PostUS until licking stops for 1 s
        SendStateMatrix(sma);
        RawEvents = RunStateMatrix;     
        
        % Update trial data -----------------------------------------------
        if ~isempty(fieldnames(RawEvents)) % If trial data was returned
            BpodSystem.Data                             = AddTrialEvents(BpodSystem.Data,RawEvents); % Computes trial events from raw data
            BpodSystem.Data.TrialSettings(currentTrial) = S; % Adds the settings used for the current trial to the Data struct (to be saved after the trial ends)
            BpodSystem.Data.TrialTypes(currentTrial)    = TrialType(currentTrial); % Adds the trial type of the current trial to data
            Outcomes                                    = UpdateOutcomePlot(TrialType, BpodSystem.Data, UsOutcome1); % update outcome plot
            BpodSystem.Data.TrialOutcome(currentTrial)  = Outcomes(currentTrial);
            SaveBpodSessionData; % Saves the field BpodSystem.Data to the current data file
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
        WaterReceived = nansum(SafeChoice)*S.RewardAmountS;
        WaterToGive   = 1000 - WaterReceived;
        disp(['WaterReceived = ' num2str(WaterReceived) 'uL']);
        disp(['WaterToGive   = ' num2str(WaterToGive) 'uL']);
    end
end

% Training stage 2 --------------------------------------------------------
if S.TrainingStage == 2 %introducing a cue and delay
    for currentTrial = 1:MaxTrials
        
        S = BpodParameterGUI('sync', S); % Synchronize the GUI
        
        % GUI Update ------------------------------------------------------
        % Controls if parameters for Tone/ Water reward are changed: if so, it is modified accordingly
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
        if ~isequal(S.TT1_Panic,TT1_PanicButon) || ~isequal(S.TT2_Panic,TT2_PanicButon)
            S.TT1_Panic = TT1_PanicButon;
            S.TT2_Panic = TT2_PanicButon;
        end
        
        % Sync parameters with BpodParameterGUI plugin
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
            if     TT1_PanicButon == 1 || Correction == 1
                TrialType(currentTrial) = 1;
            elseif TT2_PanicButon == 1 || Correction == 2
                TrialType(currentTrial) = 2;
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
        
        % Inter-trial interval --------------------------------------------
        S.ITI = 10;
        while S.ITI > 3   % ITI dustribution: 1 + exponential, truncated at 4 s
            S.ITI = exprnd(1)+1;
        end
        
        % Assemble state matrix -------------------------------------------
        sma = NewStateMatrix();
        sma = SetGlobalTimer(sma, 1, S.SoundDuration + Delay(currentTrial));
        
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
            'OutputActions', {});  % 1 + exponential foreperiod < 4s
        
        sma = AddState(sma, 'Name', 'StartStimulus', ...
            'Timer', S.SoundDuration,...
            'StateChangeConditions', {PortCondition,'Reward','Tup','Delay'},...
            'OutputActions', {'TeensyAudio1', Audio,'BNC1', 1, 'GlobalTimerTrig',1});   % play tone
        
        sma = AddState(sma, 'Name','Delay', ...
            'Timer', Delay(currentTrial),...
            'StateChangeConditions', {PortCondition,'Reward','Tup','PostUS'},...
            'OutputActions', {'PWM1', 255});
        
        sma = AddState(sma,'Name', 'Reward', ...
            'Timer',RewardValveTime,...
            'StateChangeConditions', {'Tup', 'PostUS'},...
            'OutputActions', {'ValveState', RewardValveCode, 'BNC1', 0});   % deliver water
        
        sma = AddState(sma,'Name','PostUS',...
            'Timer',1,...
            'StateChangeConditions',{PortCondition,'ResetDrinkingTimer', 'Tup','exit'},...
            'OutputActions',{'PWM1', 255});   % drinking
        
        sma = AddState(sma,'Name','ResetDrinkingTimer',...
            'Timer',0,...
            'StateChangeConditions',{'Tup','PostUS'},...
            'OutputActions',{});   % keep the animal in PostUS until licking stops for 1 s
        SendStateMatrix(sma);
        RawEvents = RunStateMatrix;
        
        % Update trial data -----------------------------------------------
        if ~isempty(fieldnames(RawEvents)) % If trial data was returned
            BpodSystem.Data                             = AddTrialEvents(BpodSystem.Data,RawEvents); % Computes trial events from raw data
            BpodSystem.Data.TrialSettings(currentTrial) = S; % Adds the settings used for the current trial to the Data struct (to be saved after the trial ends)
            BpodSystem.Data.TrialTypes(currentTrial)    = TrialType(currentTrial); % Adds the trial type of the current trial to data
            Outcomes                                    = UpdateOutcomePlot(TrialType, BpodSystem.Data, UsOutcome1); % update outcome plot
            BpodSystem.Data.TrialOutcome(currentTrial)  = Outcomes(currentTrial);
            SaveBpodSessionData; % Saves the field BpodSystem.Data to the current data file    
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
        WaterReceived = nansum(SafeChoice)*S.RewardAmountS;
        WaterToGive   = 1000 - WaterReceived;
        disp(['Water Received = ' num2str(WaterReceived) 'uL']);
        disp(['Water To Give  = ' num2str(WaterToGive) 'uL']);
    end    
end

% Training stage 3 --------------------------------------------------------
if S.TrainingStage == 3 %introducing a cue and delay
    for currentTrial = 1:MaxTrials
        
        S = BpodParameterGUI('sync', S); % Synchronize the GUI
        
        % GUI Update ------------------------------------------------------
        % Controls if parameters for Tone/ Water reward are changed: if so, it is modified accordingly
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
        if ~isequal(S.TT1_Panic,TT1_PanicButon) || ~isequal(S.TT2_Panic,TT2_PanicButon)
            S.TT1_Panic = TT1_PanicButon;
            S.TT2_Panic = TT2_PanicButon;
        end
        
        % Sync parameters with BpodParameterGUI plugin
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
            if     TT1_PanicButon == 1 || Correction == 1
                TrialType(currentTrial) = 1;
            elseif TT2_PanicButon == 1 || Correction == 2
                TrialType(currentTrial) = 2;
            else 
                TrialType(currentTrial) = ceil(rand(1,1)*2);
            end
        else 
        end
        
        % UsOutcome decision ----------------------------------------------
        clearvars A B
        if currentTrial > 20
            p = rand(1,1);
            
            A(p <= S.SafeProb) = 1;            % Safe side
            B(p <= S.RiskProb) = 2;            % Risk side - 60% big reward
            B(p > S.PunishProb & p <= 1) = 3;  % 30% punishment
            
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
            'StateChangeConditions', {'Tup', 'ITI','Port1In','RestartNoLick'},...
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
            'StateChangeConditions', {SafePort,SafePortArgument,RiskyPort,RiskyPortArgument,'Tup','Delay'},...
            'OutputActions', {'TeensyAudio1', Audio,'BNC1', 2, 'GlobalTimerTrig',1});   % play tone
        
        sma = AddState(sma, 'Name','Delay', ...
            'Timer', Delay(currentTrial),...
            'StateChangeConditions', {SafePort,SafePortArgument, RiskyPort,RiskyPortArgument,'Tup','PostUS'},...
            'OutputActions', {'BNCState', 0});
        
        sma = AddState(sma,'Name', 'RewardSafe', ...
            'Timer',RewardValveTimeSafe,...
            'StateChangeConditions', {'Tup', 'PostUS'},...
            'OutputActions', {'ValveState', RewardValveCodeSafe});   % deliver water
        
        sma = AddState(sma,'Name', 'RewardRisk', ...
            'Timer',RewardValveTimeRisk,...
            'StateChangeConditions', {'Tup', 'PostUS'},...
            'OutputActions', {'ValveState', RewardValveCodeRisk});   % deliver water
        
        sma = AddState(sma, 'Name', 'Punish', ...
            'Timer',S.PunishTime, ...
            'StateChangeConditions', {'Tup', 'PostUS'}, ...
            'OutputActions', {'ValveState', S.PunishValve}); %deiliver air puff
        
        sma = AddState(sma,'Name','PostUS',...
            'Timer',1,...
            'StateChangeConditions',{SafePort,'ResetDrinkingTimer', RiskyPort,'ResetDrinkingTimer','Tup','exit'},...
            'OutputActions',{});   % drinking
        
        sma = AddState(sma,'Name','ResetDrinkingTimer',...
            'Timer',0,...
            'StateChangeConditions',{'Tup','PostUS'},...
            'OutputActions',{});   % keep the animal in PostUS until licking stops for 1 s
        SendStateMatrix(sma);
        RawEvents = RunStateMatrix;
        
        % Update trial data -----------------------------------------------
        if ~isempty(fieldnames(RawEvents)) % If trial data was returned
            BpodSystem.Data                             = AddTrialEvents(BpodSystem.Data,RawEvents); % Computes trial events from raw data
            BpodSystem.Data.TrialSettings(currentTrial) = S; % Adds the settings used for the current trial to the Data struct (to be saved after the trial ends)
            BpodSystem.Data.TrialTypes(currentTrial)    = TrialType(currentTrial); % Adds the trial type of the current trial to data
            Outcomes                                    = UpdateOutcomePlot2(TrialType, BpodSystem.Data, UsOutcome1, UsOutcome2); % update outcome plot
            BpodSystem.Data.TrialOutcome(currentTrial)  = Outcomes(currentTrial);
            SaveBpodSessionData; % Saves the field BpodSystem.Data to the current data file
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
        WaterReceived = nansum(RiskyChoice)*S.RewardAmountB + nansum(SafeChoice)*S.RewardAmountS;
        WaterToGive   = 1000 - WaterReceived;
        disp(['Water Received = ' num2str(WaterReceived) 'uL']);
        disp(['Water To Give  = ' num2str(WaterToGive) 'uL']);
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
OutcomePlot_Gambling12(BpodSystem.GUIHandles.OutcomePlot,'update',Data.nTrials+1,1-TrialType,Outcomes, UsOutcome)
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
OutcomePlot_Gambling(BpodSystem.GUIHandles.OutcomePlot,'update',Data.nTrials+1,1-TrialType, 1-Choice ,Outcomes, [UsOutcome1;UsOutcome2])
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
