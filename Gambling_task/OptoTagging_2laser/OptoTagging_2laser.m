function OptoTagging_2laser

% BPOD protocol to perform laser tagging using 2 lasers controled by
% PulsePal, currently BLUE and YELLOW. Each laser is connected to a
% different trigger channel. Check TAGGINGBLUE and TAGGINGYELLOW for more
% information about the parameter matrix for PulsePal. 
%
% Upon starting the protocol you can choose to use ONE (B o Y) or both
% lasers. You can also choose the number of rounds to run the protocol (1
% to 5 times). If you choose both, first choose one (B or Y) and after it
% finishes choose the other (Y or B).
%
% To start the tagging click on PortIn 4. The tagging can be interrupted if
% this port is turned on again.
%
% -------------------------------------------------------------------------
% Cecília Pardo-Bellver, 2021
% Laboratory of Network Neurophysiology
% Institute of Experimental Medicine, Hungary.
% -------------------------------------------------------------------------

% Load Bpod variables
global BpodSystem

% Define parameters -------------------------------------------------------
S = BpodSystem.ProtocolSettings; %#ok<*NASGU> % Load settings chosen in launch manager into current workspace as a struct called S
S = struct;

% Read PulsePal -----------------------------------------------------------
PulsePalCOM = FindPulsePal;
PulsePal(PulsePalCOM{1,1});

% Define Lasers to do the tagging, load parameters ------------------------
Tagging = questdlg('Choose Tagging Protocol',...
    '',...
    'Blue', 'Yellow', 'Both', 'Both');

switch Tagging
    case 'Blue'
        NumberOfTags = 1;
        Tagging      = 1;
        load('TaggingBlue.mat','ParameterMatrix'); %Your ParameterMatrix for the 4 channels
        ProgramPulsePal(ParameterMatrix)
    case 'Yellow'
        NumberOfTags = 1;
        Tagging      = 2;
        load('TaggingYellow.mat','ParameterMatrix');
        ProgramPulsePal(ParameterMatrix)
    case 'Both'
        NumberOfTags = 2;
        Tagging      = 3;
end
list = {' 1 ',' 2 ',' 3 ',' 4 ',' 5 '};
MaxTrials = listdlg('PromptString','Select Number of Trials',...
    'SelectionMode','single',...
    'ListString',list);

for NumTags = 1:NumberOfTags
    if     Tagging == 1
        global PulsePalSystem
        ProgramPulsePalParam(1, 12, 1); % Channel 1(B), Responds Trigger1, Back to 1
    elseif Tagging == 2
        global PulsePalSystem %#ok<REDEF,*TLEV>
        ProgramPulsePalParam(2, 12, 1); % Channel 2(Y), Responds Trigger1, Back to 1
    elseif Tagging == 3
        laserID = questdlg('Choose Laser',...
                '',...
                'Blue', 'Yellow', 'Yellow');
        switch laserID
            case 'Blue'
                load('TaggingBlue.mat','ParameterMatrix'); %Your ParameterMatrix for the 4 channels
                ProgramPulsePal(ParameterMatrix)
                global PulsePalSystem
                ProgramPulsePalParam(1, 12, 1); % Channel 1(B), Responds Trigger1, Back to 1
            case 'Yellow'
                load('TaggingYellow.mat','ParameterMatrix');
                ProgramPulsePal(ParameterMatrix)
                global PulsePalSystem
                ProgramPulsePalParam(2, 12, 1); % Channel 2(Y), Responds Trigger1, Back to 1
        end
    end
    
    for currentTrial = 1:MaxTrials
        sma = NewStateMatrix(); % Assemble state matrix
        sma = AddState(sma,'Name', 'Start', ...
            'Timer', 1,...
            'StateChangeConditions', {'Port4In', 'PulsePalTagging'},...
            'OutputActions', {'PWM4', 255}); %Light On
        sma = AddState(sma,'Name', 'PulsePalTagging', ...
            'Timer', 120,...
            'StateChangeConditions', {'Port4In', 'Return', 'Tup', 'Return'},...
            'OutputActions', {'BNC1', 1});
        sma = AddState(sma,'Name', 'Return', ...
            'Timer', 0,...
            'StateChangeConditions', {'Tup', 'exit'},...
            'OutputActions', {'BNC1', 0});
        SendStateMatrix(sma);
        
        RawEvents = RunStateMatrix;
        if ~isempty(fieldnames(RawEvents)) % If trial data was returned
            BpodSystem.Data = AddTrialEvents(BpodSystem.Data,RawEvents); % Computes trial events from raw data
            BpodSystem.Data.TrialSettings(currentTrial) = S; % Adds the settings used for the current trial to the Data struct (to be saved after the trial ends)
            SaveBpodSessionData; % Saves the field BpodSystem.Data to the current data file
            
        end
        HandlePauseCondition;
        
    end
end
        if BpodSystem.Status.BeingUsed == 0 % If protocol was stopped, exit the loop
            EndPulsePal;
            return
        elseif currentTrial == MaxTrials
            EndPulsePal;
            return
        end
end