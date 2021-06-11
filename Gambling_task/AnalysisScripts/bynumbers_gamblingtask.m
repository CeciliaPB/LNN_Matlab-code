function StageNumbers = bynumbers_gamblingtask(animalName,sessionID,saveData)

%   Designed as an offline analysis tool for
%   behavior recorded through Bpod system, which can be executed in an unsupervised
%   manner on a daily bases. It gives a quick overview of the experiment.
%   
%   bynumbers_gamblingtask(animalName,sessionID) performs the
%   analysis for a session specified by the first two input arguments.
%
% -------------------------------------------------------------------------
% Cecília Pardo-Bellver, 2021
% Laboratory of Network Neurophysiology
% Institute of Experimental Medicine, Hungary.
% -------------------------------------------------------------------------

% Animal, session, Behavior
if nargin < 1
    disp('ERROR: No AnimalName provided.');
else
    if ischar(animalName)
        animalID = animalName;
    else
        animalID = ['BAd' num2str(animalName,'%02.0f')];
    end
end

if nargin < 2
    disp('ERROR: No sessionID provided.');
end

if nargin < 3
    saveData = 0;
end

fullpth = fullfile(getpref('cellbase','datapath'),animalID,sessionID);
cd(fullpth)

% Create Trial Events structure -------------------------------------------
if exist('TE_behaviour.mat','file') ~= 0
    TE_behaviour = load('TE_behaviour.mat');
elseif exist('TE_behaviour.mat','file') == 0
    behfile = dir([fullpth filesep animalID '_GamblingTask_'...
        sessionID(1:end-1) '*.mat']);
    behfile = behfile.name;
    TE_behaviour = TE_gamblingtask(fullfile(fullpth,behfile),1);
end

% Behavior ----------------------------------------------------------------
switch (TE_behaviour.TrainingStage(1,1))
    case 1 
        % Calculate the mean time the animals need to lick to each side.
        % Criterion for next stage: time < Ns in seconds in both sides.
        
        Ns = 5;
        TrialType  = TE_behaviour.TrialType; % 1 = R/ LargeR; 2 = L/ LargeL
        
        % Licks R
        NumTrials = length(TE_behaviour.LickRIn);
        LickRIn = TE_behaviour.LickRIn;
        Rtimes = nan(1,NumTrials-1);
        for ii = 5:NumTrials-1 % The last trial might not be completed
            if isempty(LickRIn{ii})
            elseif  TrialType(ii)==1
                Rtimes(ii) = LickRIn{ii}(1);
            end
        end
        mdTimeR = nanmedian(Rtimes);
        
        % Licks L
        NumTrials = length(TE_behaviour.LickLIn);
        LickLIn = TE_behaviour.LickLIn;
        Ltimes = nan(1,NumTrials-1);
        for ii = 5:NumTrials-1 % The last trial might not be completed
            if isempty(LickLIn{ii})
            elseif TrialType(ii)==2
                Ltimes(ii) = LickLIn{ii}(1);
            end
        end
        mdTimeL = nanmedian(Ltimes);
        
        if mdTimeR<Ns && mdTimeL<Ns
            NextStg = 1;
        else 
            NextStg = 0;
        end
        
        StageNumbers.Stage   = TE_behaviour.TrainingStage(1,1);
        StageNumbers.mdTimeR = mdTimeR;
        StageNumbers.mdTimeL = mdTimeL;
        StageNumbers.NextStg = NextStg;
        
    case 2
        % Calculate the percentage of correct first licks to each side.
        % Criterion for next stage: percentage > Correct100 in both sides.
        
        Correct100 = 50;
        TrialType  = TE_behaviour.TrialType; % 1 = R/ LargeR; 2 = L/ LargeL
        FirstLick  = TE_behaviour.FirstLick; % 1 = R; 2 = L
        
        % Licks R
        nTT1 = sum(double(TrialType == 1));
        RlickTT1 = sum(double(FirstLick == 1 & TrialType == 1));
        CorrectR = RlickTT1/nTT1*100;
        
        % Licks L
        nTT2 = sum(double(TrialType == 2));
        LlickTT2 = sum(double(FirstLick == 2 & TrialType == 2));
        CorrectL = LlickTT2/nTT2*100;
        
        if CorrectR>Correct100 && CorrectL>Correct100
            NextStg = 1;
        else 
            NextStg = 0;
        end
        
        StageNumbers.Stage       = TE_behaviour.TrainingStage(1,1);
        StageNumbers.CorrectR100 = CorrectR;
        StageNumbers.CorrectL100 = CorrectL;
        StageNumbers.NextStg     = NextStg;
        
    case 3
        NumTrials  = length(TE_behaviour.TrialType);
        TrialType  = TE_behaviour.TrialType; % 1 = R/ LargeR; 2 = L/ LargeL
        SafeChoice = TE_behaviour.SafeChoice;
        RiskChoice = TE_behaviour.RiskyChoice;
        Punishment = TE_behaviour.Punishment;
        FirstLick  = TE_behaviour.FirstLick; % 1 = R; 2 = L
        nTT2 = sum(double(TrialType == 2));
        nTT1 = sum(double(TrialType == 1));
        
        % Licks R
        LlickTT1 = sum(double(SafeChoice == 1 & TrialType == 1));
        RlickTT1 = sum(double(RiskChoice == 1 & TrialType == 1));
        Rpunish  = sum(double(Punishment == 1));
        Rlick    = (sum(double(FirstLick == 1)) + Rpunish)/NumTrials*100;
        Rrisk    = RlickTT1/nTT1*100;
        Rsafe    = LlickTT1/nTT1*100;
        
        % Licks L
        LlickTT2 = sum(double(RiskChoice == 1 & TrialType == 2));
        RlickTT2 = sum(double(SafeChoice == 1 & TrialType == 2));
        Lpunish  = sum(double(Punishment == 2));
        Llick    = (sum(double(FirstLick == 2))+ Lpunish)/NumTrials*100;
        Lsafe    = RlickTT2/nTT2*100;
        Lrisk    = LlickTT2/nTT2*100;
        
        StageNumbers.Stage    = TE_behaviour.TrainingStage(1,1);
        StageNumbers.Rrisk100 = Rrisk;
        StageNumbers.Rsafe100 = Rsafe;
        StageNumbers.Lrisk100 = Lrisk;
        StageNumbers.Lsafe100 = Lsafe;
        StageNumbers.Punish100   = nansum(Punishment)/(nTT1 + nTT2)*100;
        StageNumbers.Rpunish100  = Rpunish/nTT1*100;
        StageNumbers.Lpunish100  = Lpunish/nTT2*100;
        StageNumbers.Rlick100 = Rlick;
        StageNumbers.Llick100 = Llick;
end

        if saveData ~= 0
            savename = 'StageByNumbers.mat';
            save([fullpth filesep savename],'-struct','StageNumbers')
        end
        
end