function Poking2AFC_3_Test

%  Freely inspired by the Jaramillo & Zador 2014 task (doi:
% 10.3389/fnsys.2014.00173), making it a working memory task and not an
% attention one. Substantially improved thanks to Juliane Martin from Duda lab and Sümegi Mate from the Nusser one.
%
% Third step: are the mice ready to rocksteady? Central reward is now composed of 3 blocks
% (S.GUI.FixationTime x 3, with micro reward for engegement in the centre), while you
% decide the total length of the CS.
% Basically what you do now is reducing day by day the CS length from
% infinity to 0.2, at which point it will be shorter than the mandatory
% fixation time length.
% There is airpuff.
%
%
% Nicola Solari, 2018
% Lendulet Laboratory of Systems Neuroscience (hangyalab.koki.hu)
% Institute of Experimental Medicine, Hungarian Academy of Sceinces

clc
global BpodSystem

prompt = {'Cue_Length'};
title = 'CS Shortening';
dims = [1 35];
answer = inputdlg(prompt,title,dims); % Insert the value of the CS that you want the mouse to have
Cue_Length = str2num(answer{1});

% Structure creation
S = BpodSystem.ProtocolSettings; % Load settings chosen in launch manager into current workspace as a struct called S
S = struct;
if isempty(fieldnames(S))  % If settings file was an empty struct, populate struct with default settings
    S.GUI.NumTrialTypes = 2;   % Different trial types corresponding to different cues and outcome contingencies
    S.GUI.FixationTime = 0.1; % THIS IS MULTIPLIED BY 3!!!!
    S.GUI.TrialWindow = 4; % Once central reward has been delivered, the mouse has a fixed time window to make his choice
    S.GUI.SinWavekHz1 = 3; % Cue tone #1 in kHz - tone #1 signals the availability of the reward in the left arm
    S.GUI.SinWavedB1 = 72; % Cue tone #1 dB SPL
    S.GUI.SinWavekHz2 = 11; % Cue tone #2 in kHz - tone #2 signals the availability of the reward in the right arm
    S.GUI.SinWavedB2 = 72; % Cue tone #2 dB SPL
    S.GUI.WhiteWavedB = 95; % WhiteNoise tone dB SPL
    S.GUI.RewardObtained = 0;
    S.GUI.nITI = 1; % Normal Inter Trial Interval
    S.GUI.pITI = 10; % Punishment Inter Trial Interval
    S.GUI.LeftPanicButon = 0; % For Trial Type 1 Imposition
    S.GUI.RightPanicButon = 0; % For Trial Type 2 Imposition
    S.GUI.Type1_Success = 0; % Live monitoring of performance in Trial Type 1
    S.GUI.Type2_Success = 0; % Live monitoring of performance in Trial Type 2
    S.GUI.Water_Left = 2; % Reward to the left port (it's interactive, you can change it during the task)
    S.GUI.Water_Right = 2; % Reward to the right port (it's interactive, you can change it during the task)
end
BpodParameterGUI('init', S); % Initialize parameter GUI plugin
LPB = S.GUI.LeftPanicButon;
RPB = S.GUI.RightPanicButon;
WatL = S.GUI.Water_Left;
WatR = S.GUI.Water_Right;

% Trials definition
MaxTrials = 800; % Maximum number of trials
rng('shuffle')   % Reset pseudorandom seed
BpodSystem.Data.TrialTypes = []; % BpodSystem.Data.TrialTypes = []
TrialTypes = [];

% Teensy & calibration data engagement
T = TeensyAudioPlayer('COM1'); % Opens Teensy Server at the COM port engaged by the device
FilePath = fullfile(BpodSystem.Path.LocalDir,'Protocols','TeensyCalibration','TeensyCalData.mat');
load(FilePath); % Load Calibration data as a reference

% Tones 1 and 2 Creation
Tg1 = S.GUI.SinWavedB1; % Wanted dB for tone 1
Tg2 = S.GUI.SinWavedB2; % Wanted dB for tone 2
Tg99 =  S.GUI.WhiteWavedB;
Fr1 = S.GUI.SinWavekHz1;
Fr2 = S.GUI.SinWavekHz2;
Fr99 = 22; % Becouse the WhiteNoise is saved at the 22nd during the calibration (for me, having done it for 21 tracks)
SPL1 = TeensyCalData.SPL(Fr1); % Recalls calibrated dB for the frequency of tone 1
SPL2 = TeensyCalData.SPL(Fr2); % Recalls calibrated dB for the frequency of tone 2
SPL99 = TeensyCalData.SPL(Fr99); % Recalls calibrated dB for the WhiteNoise
Ampl1 = TeensyCalData.Amplitude(Fr1); % Recalls calibrated amplitude for the tone 1 frequency
Ampl2 = TeensyCalData.Amplitude(Fr2); % Recalls calibrated amplitude for the tone 2 frequency
Ampl99 = TeensyCalData.Amplitude(Fr99); % Recalls calibrated amplitude for the WhiteNoise
NewAmpl1  = AmplAdjst(SPL1,Tg1,Ampl1); % Calculates new amplitude for tone 1
NewAmpl2  = AmplAdjst(SPL2,Tg2,Ampl2); % Calculates new amplitude for tone 2
NewAmpl99  = AmplAdjst(SPL99,Tg99,Ampl99); % Calculates new amplitude for WhiteNoise
sinewave1  = NewAmpl1.*sin(2*pi*Fr1*1000/44100.*(0:44100*Cue_Length)); % Creates the long sinewave of tone 1
sinewave2  = NewAmpl2.*sin(2*pi*Fr2*1000/44100.*(0:44100*Cue_Length)); % Creates the long sinewaves of tone 2
sinewave99  = NewAmpl99*randn(1,2*44100); % Creates the long sinewaves of WhiteNoise
T.load(Fr1, sinewave1); % Upload to SD card the sinwave
T.load(Fr2, sinewave2); % Upload to SD card the sinwave
T.load(99, sinewave99); % Upload to SD card the sinwave

CentreReward = 0.2; % Central valve micro reward amount
R3 = GetValveTimes(CentreReward, 3);
CentreValveTime = R3;

for currentTrial= 1:800
    % Controls if a Panic Button has been activated: if so, the corresponding Trial Type is imposed untill counterorder
    cLPB = S.GUI.LeftPanicButon;
    cRPB = S.GUI.RightPanicButon;
    if ~isequal(LPB,cLPB) || ~isequal(RPB,cRPB)
        if cLPB == 0
            LeftSwitch = 0;
        else
            LeftSwitch = 1;
        end
        if cRPB == 0
            RightSwitch = 0;
        else
            RightSwitch = 1;
        end
        LPB = cLPB;
        RPB = cRPB;
    end
    
    % Check if water reward has changed, and if so adjust valve delivery time
    cWatL = S.GUI.Water_Left;
    cWatR = S.GUI.Water_Right;
    if ~isequal(WatL,cWatL) || ~isequal(WatR,cWatR)
        WatL = cWatL;
        WatR = cWatR;
    end
    LeftReward = WatL;
    RightReward = WatR;
    R1 = GetValveTimes(WatL, 1); LeftValveTime = R1;
    R2 = GetValveTimes(WatR, 3);  RightValveTime = R2;
    
    % Synchronize the GUI!!!! Put it or it won't work!!!
    S = BpodParameterGUI('sync', S);
    
    % Trial Type decision
    if currentTrial<=16 % first 16 trials are purely random
        Correction=0;
        LeftSwitch = 0;
        RightSwitch = 0;
    else
        Correction=BiasCounter(TrialTypes, BpodSystem.Data);
    end
    if    LeftSwitch == 1 || Correction == 1
        TrialTypes(currentTrial) = 1;
        C_Tag = 'Corr1';
    elseif   RightSwitch == 1 || Correction == 2
        TrialTypes(currentTrial) = 2;
        C_Tag = 'Corr2';
    else
        TrialTypes(currentTrial) = ceil(rand(1,1)*2);
        C_Tag = 'Corr0';
    end
    disp(['Trial# ' num2str(currentTrial) ' TrialType: ' num2str(TrialTypes(currentTrial))]);
    
    % Performance and Step Selection
    if currentTrial == 1;
        Performance = 0;
    end
    
    % Trial Types Parameters
    if TrialTypes(currentTrial) == 1
        CorrectPort = 'Port1In'; UncorrectPort = 'Port2In';
        Reward = 'LMajorReward'; Punishment = 'LFail';
        MajorValve = 1; MajorTime = LeftValveTime; Audio = Fr1;
    else
        CorrectPort = 'Port2In'; UncorrectPort = 'Port1In';
        Reward = 'RMajorReward'; Punishment = 'RFail';
        MajorValve = 2; MajorTime = RightValveTime; Audio = Fr2;
    end
    
    Performance % leave uncommented to see how many good trials the mouse did
    
    % Assemble state matrix
    sma = NewStateMatrix();
    sma = SetGlobalTimer(sma, 'TimerID', 1, 'Duration', S.GUI.TrialWindow);
    % the tags (for fine-grain analysis to distinguish "natural" and "imposed" Trial Types)
    sma = AddState(sma, 'Name', 'C_Tag', ... % Find Correction value and select tag event
        'Timer',0,...
        'StateChangeConditions', {'Tup',C_Tag},...
        'OutputActions',  {});
    sma = AddState(sma, 'Name', 'Corr1', ... % Event Tag for Correction 1
        'Timer',0,...
        'StateChangeConditions', {'Tup','WaitForPoking0'},...
        'OutputActions',  {});
    sma = AddState(sma, 'Name', 'Corr2', ... % Event Tag for Correction 2
        'Timer',0,...
        'StateChangeConditions', {'Tup','WaitForPoking0'},...
        'OutputActions',  {});
    sma = AddState(sma, 'Name', 'Corr0', ... % Event Tag for Correction 0
        'Timer',0,...
        'StateChangeConditions', {'Tup','WaitForPoking0'},...
        'OutputActions',  {});
    % the stem
    sma = AddState(sma, 'Name', 'WaitForPoking0', ...  % For Bpod/OE syncronization
        'Timer', 0.01,...
        'StateChangeConditions', { 'Tup', 'WaitForPoking'},...
        'OutputActions', {'BNCState', 1});
    sma = AddState(sma, 'Name', 'WaitForPoking', ...  % Mouse has to poke centre to proceed
        'Timer', 0,...
        'StateChangeConditions', { 'Port3In', 'PokeInLight_1'},...
        'OutputActions', {'LED', 3,'PWM8', 80});
    sma = AddState(sma, 'Name', 'PokeInLight_1', ... % Light in the central port and cue
        'Timer', S.GUI.FixationTime,...
        'StateChangeConditions', {'Tup','PokeInLight_2','Port3Out','aITI'},...
        'OutputActions', {'LED', 3,'TeensyAudio1', Audio,'PWM8', 80});
    sma = AddState(sma, 'Name', 'PokeInLight_2', ... % Light in the central port, 1st water
        'Timer', CentreValveTime,...
        'StateChangeConditions', {'Tup','PokeInLight_3','Port3Out','aaITI'},...
        'OutputActions', {'LED', 3,'ValveState', 4,'PWM8', 80});
    sma = AddState(sma, 'Name', 'PokeInLight_3', ... % Light in the central port
        'Timer', S.GUI.FixationTime,...
        'StateChangeConditions', {'Tup','PokeInLight_4','Port3Out','aITI'},...
        'OutputActions', {'LED', 3,'PWM8', 80});
    sma = AddState(sma, 'Name', 'PokeInLight_4', ... % Light in the central port, 2nd water
        'Timer', CentreValveTime,...
        'StateChangeConditions', {'Tup','PokeInLight_5','Port3Out','aaITI'},...
        'OutputActions', {'LED', 3,'ValveState', 4,'PWM8', 80});
    sma = AddState(sma, 'Name', 'PokeInLight_5', ... % Light in the central port
        'Timer', S.GUI.FixationTime,...
        'StateChangeConditions', {'Tup','PokeInLight_6','Port3Out','aITI'},...
        'OutputActions', {'LED', 3,'PWM8', 80});
    sma = AddState(sma, 'Name', 'PokeInLight_6', ... % Light in the central port, 3rd water
        'Timer', CentreValveTime,...
        'StateChangeConditions', {'Tup','TempusFit','Port3Out','aaITI'},...
        'OutputActions', {'LED', 3,'ValveState', 4,'PWM8', 80});
    sma = AddState(sma, 'Name', 'TempusFit',... % Start Timer
        'Timer', 0,...
        'StateChangeConditions', {'Tup','LuxFit'},...
        'OutputActions', {'GlobalTimerTrig',1,'PWM8', 80});
    % the choice
    sma = AddState(sma, 'Name', 'LuxFit', ... % Light on the sides
        'Timer', 0,...
        'StateChangeConditions', {CorrectPort, Reward, UncorrectPort, Punishment,'GlobalTimer1_End','tITI'},...
        'OutputActions', {'LED', 1,'LED', 2,'PWM8', 80});
    sma = AddState(sma, 'Name', 'LMajorReward', ... % Mouse get the major reward on the Left
        'Timer', LeftValveTime,...
        'StateChangeConditions', {'Tup', 'ITI'},...
        'OutputActions', {'TeensyAudio1', 254,'ValveState', 1,'PWM8', 80});
    sma = AddState(sma, 'Name', 'RMajorReward', ... % Mouse get the major reward on the Right
        'Timer', RightValveTime,...
        'StateChangeConditions', {'Tup', 'ITI'},...
        'OutputActions', {'TeensyAudio1', 254,'ValveState', 2,'PWM8', 80});
    sma = AddState(sma, 'Name', 'ITI', ... % Breath and restart all over agin
        'Timer', S.GUI.nITI,...
        'StateChangeConditions', {'Tup', 'exit'},...
        'OutputActions', {'PWM8', 80});
    sma = AddState(sma, 'Name', 'LFail', ... % Mouse fails on the Left
        'Timer', 0.3,...
        'StateChangeConditions', {'Tup', 'pITI'},...
        'OutputActions', {'ValveState',128,'TeensyAudio1',99,'PWM8', 80}); %
    sma = AddState(sma, 'Name', 'RFail', ... % Mouse fails on the Right
        'Timer', 0.3,...
        'StateChangeConditions', {'Tup', 'pITI'},...
        'OutputActions', {'ValveState',128,'TeensyAudio1',99,'PWM8', 80}); %
    sma = AddState(sma, 'Name', 'pITI', ... % The ITI of the wrong choice aftert the central poking
        'Timer', S.GUI.pITI,...
        'StateChangeConditions', {'Tup', 'exit'},...
        'OutputActions', {'PWM8', 80});
    sma = AddState(sma, 'Name', 'tITI', ... % Time Out: mouse did not choose
        'Timer', S.GUI.pITI,...
        'StateChangeConditions', {'Tup', 'exit'},...
        'OutputActions', {'TeensyAudio1', 254,'PWM8', 80});
    sma = AddState(sma, 'Name', 'aITI', ... % Abort: mice withdrawed during the mandatory fixation period (light only)
        'Timer', S.GUI.pITI,...
        'StateChangeConditions', {'Tup', 'exit'},...
        'OutputActions', {'TeensyAudio1', 99,'PWM8', 80});
    sma = AddState(sma, 'Name', 'aaITI', ... % Abort: mice withdrawed during the mandatory fixation period (while micro-reward was delivered)(yes, is a suprefluous matrix)
        'Timer', S.GUI.pITI,...
        'StateChangeConditions', {'Tup', 'exit'},...
        'OutputActions', {'TeensyAudio1', 99, 'PWM8', 80});
    SendStateMatrix(sma);
    RawEvents = RunStateMatrix;
    if ~isempty(fieldnames(RawEvents)) % If trial data was returned
        BpodSystem.Data = AddTrialEvents(BpodSystem.Data,RawEvents); % Computes trial events from raw data
        BpodSystem.Data.TrialSettings(currentTrial) = S; % Adds the settings used for the current trial to the Data struct (to be saved after the trial ends)
        BpodSystem.Data.TrialTypes(currentTrial) = TrialTypes(currentTrial); % Adds the trial type of the current trial to data
        SaveBpodSessionData; % Saves the field BpodSystem.Data to the current data file
    end
    UpdateWaterReward(currentTrial, BpodSystem.Data, RightReward, LeftReward, CentreReward);
    [L_performance, R_performance] = UpdateSidePerformance(currentTrial, BpodSystem.Data);
    Performance = CheckSuccess(TrialTypes, BpodSystem.Data);
    HandlePauseCondition;
    if BpodSystem.Status.BeingUsed == 0
        return
    end
end

function [CalAmpl] = AmplAdjst(SPL,Tg,Ampl) % Calculate the new proper sinewave amplitude
y = SPL - Tg;
b =  20 * log10(Ampl) - y;
c = b / 20;
CalAmpl = 10 .^ c;

function  UpdateWaterReward(currentTrial, Data, RightReward, LeftReward, CentreReward) % Sums up the water obtained    [TotalReward] =
global BpodSystem
Outcomes = zeros(1,Data.nTrials);
for x = 1:Data.nTrials
    if ~isnan(Data.RawEvents.Trial{x}.States.LMajorReward(1))
        Outcomes(x) = LeftReward + CentreReward;
    elseif ~isnan(Data.RawEvents.Trial{x}.States.RMajorReward(1))
        Outcomes(x) = RightReward + CentreReward;
    elseif ~isnan(Data.RawEvents.Trial{x}.States.RFail(1))
        Outcomes(x) = CentreReward;
    elseif ~isnan(Data.RawEvents.Trial{x}.States.LFail(1))
        Outcomes(x) = CentreReward;
    else
        Outcomes(x) = 0;
    end
    TotalReward = sum(Outcomes);
end
NumParam = BpodSystem.GUIData.ParameterGUI.nParams;
for  n = 1:NumParam;
    if strcmp(BpodSystem.GUIData.ParameterGUI.ParamNames(n), 'RewardObtained') == 1
        b = n;
    end
end
if currentTrial ~= 1
    ThisParamHandle = BpodSystem.GUIHandles.ParameterGUI.Params(b);
    set(ThisParamHandle, 'String', num2str(TotalReward));
end

function [Correction] = BiasCounter(TrialTypes, Data) % mild correction for bias (...maybe.)
global BpodSystem
window = 15;
% Outcomes = zeros(1,1000);
for x = 1:Data.nTrials
    if ~isnan(Data.RawEvents.Trial{x}.States.LFail(1)) % wrong TrialType1 trial
        Outcomes(x) = 1;
    elseif ~isnan(Data.RawEvents.Trial{x}.States.RFail(1)) % wrong TrialType2 trial
        Outcomes(x) = 2;
    elseif ~isnan(Data.RawEvents.Trial{x}.States.aITI(1)) % abort (both)
        Outcomes(x) = 9    ;
    elseif ~isnan(Data.RawEvents.Trial{x}.States.aaITI(1)) % abort (both)
        Outcomes(x) = 9;
    else
        Outcomes(x) = 3; % the rewards or the failure
    end
    if x > window
        LastOutcomes = Outcomes(end-14:end);
        t = tabulate(LastOutcomes);
        if x > 15
            if  t(1,3) >50
                C = 1
            elseif  t(2,3) >50
                C = 2
            else
                C = 0
            end
        else
        end
        Correction = C;
    end
end

function [Performance] = CheckSuccess(TrialTypes, Data) % self explanatory
global BpodSystem
Outcomes = zeros(1,Data.nTrials);
for x = 1:Data.nTrials
    if ~isnan(Data.RawEvents.Trial{x}.States.LMajorReward(1))
        Outcomes(x) = 1;
    elseif ~isnan(Data.RawEvents.Trial{x}.States.RMajorReward(1))
        Outcomes(x) = 1;
    else
        Outcomes(x) = 0;
    end
    Performance = sum(Outcomes);
end

function [L_performance, R_performance] = UpdateSidePerformance(currentTrial, Data) % Calculates the side performance
global BpodSystem
L_fail = zeros (1,700);
R_fail = zeros (1,700);
L_success = zeros (1,700);
R_success = zeros (1,700);
for x = 1:Data.nTrials
    if ~isnan(Data.RawEvents.Trial{x}.States.LFail(1)) % wrong TrialType1 trial
        Endings(x) = 1
    elseif ~isnan(Data.RawEvents.Trial{x}.States.RFail(1))
        Endings(x) = 2
    elseif ~isnan(Data.RawEvents.Trial{x}.States.LMajorReward(1))
        Endings(x) = 3
    elseif ~isnan(Data.RawEvents.Trial{x}.States.RMajorReward(1))
        Endings(x) = 4
    else
        Endings(x) = 5
    end
    L_fail = sum(Endings == 1) % the Trial Type 1 Failure
    R_fail = sum(Endings == 2) % the Trial Type 2 Failure
    L_success = sum(Endings == 3) % the Trial Type 1 Reward
    R_success = sum(Endings == 4) % the Trial Type 2 Reward
    L_performance = L_success/(L_success+L_fail)
    R_performance = R_success/(R_success+R_fail)
end
NumParam = BpodSystem.GUIData.ParameterGUI.nParams;
for  n = 1:NumParam;
    if strcmp(BpodSystem.GUIData.ParameterGUI.ParamNames(n), 'Type1_Success') == 1
        b = n;
    end
end
if currentTrial ~= 1
    ThisParamHandle = BpodSystem.GUIHandles.ParameterGUI.Params(b);
    set(ThisParamHandle, 'String', num2str(L_performance));
end
NumParam = BpodSystem.GUIData.ParameterGUI.nParams;
for  n = 1:NumParam;
    if strcmp(BpodSystem.GUIData.ParameterGUI.ParamNames(n), 'Type2_Success') == 1
        b = n;
    end
end
if currentTrial ~= 1
    ThisParamHandle = BpodSystem.GUIHandles.ParameterGUI.Params(b);
    set(ThisParamHandle, 'String', num2str(R_performance));
end

function [Correction] = BiasCounter2(TrialTypes, Data) % not in use, but let's keep it here
global BpodSystem
window = 15;
for x = 1:Data.nTrials
    if ~isnan(Data.RawEvents.Trial{x}.States.LMajorReward(1)) % Correct TrialType1 trial
        Outcomes(x) = 1;
    elseif ~isnan(Data.RawEvents.Trial{x}.States.RMajorReward(1)) % Correct TrialType2 trial
        Outcomes(x) = 2;
    elseif ~isnan(Data.RawEvents.Trial{x}.States.LFail(1)) % Uncorrect TrialType1 trial
        Outcomes(x) = 11;
    elseif ~isnan(Data.RawEvents.Trial{x}.States.RFail(1)) % Uncorrect TrialType2 trial
        Outcomes(x) = 22;
    else
        Outcomes(x) = 9; % abort and timeout
    end
    if x > window
        LastOutcomes = Outcomes(end-14:end);
        Num_Corr_Left = length(LastOutcomes(LastOutcomes==1));
        Num_Corr_Right =  length(LastOutcomes(LastOutcomes==2));
        Num_Uncorr_Left = length(LastOutcomes(LastOutcomes==11));
        Num_Uncorr_Right = length(LastOutcomes(LastOutcomes==22));
        Perc_Uncorr_Left = Num_Uncorr_Left/(Num_Corr_Left+Num_Uncorr_Left);
        Perc_Uncorr_Right = Num_Uncorr_Right/(Num_Corr_Right+Num_Uncorr_Right);
        if x > 15
            if Perc_Uncorr_Left >= 0.5 && Perc_Uncorr_Right >= 0.5
                C = 0
            elseif Perc_Uncorr_Left <= 0.5 && Perc_Uncorr_Right <= 0.5
                C = 0
            elseif Perc_Uncorr_Left >= 0.5 && Perc_Uncorr_Right < 0.5
                C = 1
            elseif Perc_Uncorr_Left < 0.5 && Perc_Uncorr_Right >= 0.5
                C = 2
            else
                C = 0
            end
        else
            
        end
        Correction = C;
    end
end