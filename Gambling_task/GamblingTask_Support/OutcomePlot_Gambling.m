function OutcomePlot_Gambling(AxesHandle, Action, varargin)
%
% -------------------------------------------------------------------------
%{
This file is part of the Bpod Project
Copyright (C) 2014 Joshua I. Sanders, Cold Spring Harbor Laboratory, NY, USA

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, version 3.

This program is distributed  WITHOUT ANY WARRANTY and without even the 
implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
%}
% -------------------------------------------------------------------------
%
% My lines: OutcomePlot_Gambling(BpodSystem.GUIHandles.OutcomePlot,'init',[zeros(1,length(Trials)) ; ones(1,length(Trials))], [UsOutcome1;UsOutcome2]);
%           OutcomePlot_Gambling(BpodSystem.GUIHandles.OutcomePlot,'update',Data.nTrials+1,1-TrialType, 1-Choice ,Outcomes, [UsOutcome1;UsOutcome2])
%
% Plug in to Plot reward side and trial outcome.
% AxesHandle = handle of axes to plot on
% Action = specific action for plot, "init" - initialize OR "update" -  update plot
%
%Example usage:
% OutcomePlot(AxesHandle,'init',TrialTypeSides)
% OutcomePlot(AxesHandle,'init',TrialTypeSides,'ntrials',90)
% OutcomePlot(AxesHandle,'update',CurrentTrial,TrialTypeSides,OutcomeRecord)
%
% varargins:
% TrialTypeSides: Vector of 0's (right) or 1's (left) to indicate reward side (0,1), or 'None' to plot trial types individually
% OutcomeRecord:  Vector of trial outcomes
%                 Simplest case: 
%                               1: correct trial (green)
%                               0: incorrect trial (red)
%                 Advanced case: 
%                               NaN: future trial (blue)
%                                -1: withdrawal (red circle)
%                                 0: incorrect choice (red dot)
%                                 1: correct choice (green dot)
%                                 2: did not choose (green circle)
% OutcomeRecord can also be empty
% Current trial: the current trial number
%
% Adapted from BControl (SidesPlotSection.m) 
% Kachi O. 2014.Mar.17
%
% -------------------------------------------------------------------------
% Modified by Cecília Pardo-Bellver, 2020
% Laboratory of Network Neurophysiology
% Institute of Experimental Medicine, Hungary.
% -------------------------------------------------------------------------

%% Code Starts Here
global nTrialsToShow %this is for convenience
% global BpodSystem

switch Action
    case 'init'
        %initialize pokes plot
        TrialTypesSide = varargin{1};
        ColorList = varargin{2};
                       
        nTrialsToShow = 20; %default number of trials to display
        
        if nargin > 4 %custom number of trials
            nTrialsToShow =varargin{4};
        end
        
        %plot in specified axes
        scatter(AxesHandle,  find(ColorList(1,1:nTrialsToShow) == 0), zeros(1,length(find(ColorList(1,1:nTrialsToShow) == 0))),'MarkerFaceColor',[1 1 1],'MarkerEdgeColor', [0 0 1]);
        hold(AxesHandle, 'on');
        scatter(AxesHandle,  find(ColorList(2,1:nTrialsToShow) == 0), repmat(-1, 1, length(find(ColorList(2,1:nTrialsToShow) == 0))),'MarkerFaceColor',[1 1 1],'MarkerEdgeColor', [0 0 1]);
        scatter(AxesHandle,  find(ColorList(1,1:nTrialsToShow) == 1), zeros(1,length(find(ColorList(1,1:nTrialsToShow) == 1))),'MarkerFaceColor',[1 1 1],'MarkerEdgeColor', [0 1 0]);
        scatter(AxesHandle,  find(ColorList(2,1:nTrialsToShow) == 2), repmat(-1, 1, length(find(ColorList(2,1:nTrialsToShow) == 2))),'MarkerFaceColor',[1 1 1],'MarkerEdgeColor', [0 0.75 0.25]);
        scatter(AxesHandle,  find(ColorList(2,1:nTrialsToShow) == 3), repmat(-1, 1, length(find(ColorList(2,1:nTrialsToShow) == 3))),'MarkerFaceColor',[1 1 1],'MarkerEdgeColor', [1 0 0]);
        scatter(AxesHandle,  find(TrialTypesSide(1,1:nTrialsToShow) == 0), ones(1, length(find(TrialTypesSide(1,1:nTrialsToShow) == 0))),'MarkerFaceColor',[1 0.3 0],'MarkerEdgeColor', [1 0.3 0]);
        scatter(AxesHandle,  find(TrialTypesSide(1,1:nTrialsToShow) == -1), ones(1, length(find(TrialTypesSide(1,1:nTrialsToShow) == -1))),'MarkerFaceColor',[0 0.3 1],'MarkerEdgeColor', [0 0.3 1]);
        set(AxesHandle,'TickDir', 'out','XLim',[0 21],'YLim', [-2, 2], 'YTick', [-1 0 1],'YTickLabel', {'Risk','Safe', 'TrialType'}, 'FontSize', 10);
        xlabel(AxesHandle, 'Trial#', 'FontSize', 18);
        
        
    case 'update'
        CurrentTrial = varargin{1};
        TrialTypesSide = varargin{2};
        Choice = varargin{3};
        OutcomeRecord = varargin{4};
        ColorList = varargin{5};
       
        if CurrentTrial<1
            CurrentTrial = 1;
        end
        
        % recompute xlim
        [mn, ~] = rescaleX(AxesHandle,CurrentTrial,nTrialsToShow);
        
        %plot future trials
        NewTrialsIndx = mn:CurrentTrial-1;
        scatter(AxesHandle,  NewTrialsIndx(ColorList(1,NewTrialsIndx) == 0), zeros(1, length(NewTrialsIndx(ColorList(1,NewTrialsIndx) == 0))),'MarkerFaceColor',[1 1 1],'MarkerEdgeColor', [0 0 1]);
        hold(AxesHandle, 'on');
        scatter(AxesHandle,  NewTrialsIndx(ColorList(2,NewTrialsIndx) == 0), repmat(-1,1,length(NewTrialsIndx(ColorList(2,NewTrialsIndx) == 0))),'MarkerFaceColor',[1 1 1],'MarkerEdgeColor', [0 0 1]);

        scatter(AxesHandle,  NewTrialsIndx(ColorList(1,NewTrialsIndx) == 1), zeros(1,length(NewTrialsIndx(ColorList(1,NewTrialsIndx) == 1))),'MarkerFaceColor',[1 1 1],'MarkerEdgeColor', [0 1 0]);
        scatter(AxesHandle,  NewTrialsIndx(ColorList(2,NewTrialsIndx) == 2), repmat(-1,1,length(NewTrialsIndx(ColorList(2,NewTrialsIndx) == 2))),'MarkerFaceColor',[1 1 1],'MarkerEdgeColor', [0 0.75 0.25]);
        scatter(AxesHandle,  NewTrialsIndx(ColorList(2,NewTrialsIndx) == 3), repmat(-1,1,length(NewTrialsIndx(ColorList(2,NewTrialsIndx) == 3))),'MarkerFaceColor',[1 1 1],'MarkerEdgeColor', [1 0 0]);
        
%         %Plot current trial
%         scatter(AxesHandle,CurrentTrial,SideList(CurrentTrial), 'o', 'MarkerFaceColor',[1 1 1],'MarkerEdgeColor', 'k')
%         scatter(AxesHandle,CurrentTrial,SideList(CurrentTrial), '+', 'MarkerFaceColor',[1 1 1],'MarkerEdgeColor', 'k')
        
        %Plot past trials
        if ~isempty(OutcomeRecord)
            indxToPlot = mn:CurrentTrial-1;
            
            % Plot TrialType
            TT1 = (TrialTypesSide(indxToPlot) == 0);
            TT2 = (TrialTypesSide(indxToPlot) == -1);
            scatter(AxesHandle,  indxToPlot(TT1), ones(1,length(indxToPlot(TT1))),'MarkerFaceColor',[1 0.3 0],'MarkerEdgeColor', [1 0.3 0]);
            scatter(AxesHandle,  indxToPlot(TT2), ones(1,length(indxToPlot(TT2))),'MarkerFaceColor',[0 0.3 1],'MarkerEdgeColor', [0 0.3 1]);
        
            %Plot NoLick, NoReward
            DidNotChooseTrialsIndx   = (OutcomeRecord(indxToPlot) == 0);
            scatter(AxesHandle,  indxToPlot(DidNotChooseTrialsIndx),...
                Choice(DidNotChooseTrialsIndx),'bo','MarkerFaceColor',[1 1 1]);
            %Plot Safe lick, rewarded
            SafeCorrectTrialsIndx    = (OutcomeRecord(indxToPlot) == 1);
            scatter(AxesHandle,  indxToPlot(SafeCorrectTrialsIndx),...
                Choice(indxToPlot(SafeCorrectTrialsIndx)),'MarkerFaceColor','g','MarkerEdgeColor', 'g');
            %Plot Risk lick, rewarded
            RiskCorrectTrialsIndx    = (OutcomeRecord(indxToPlot) == 2);
            scatter(AxesHandle,  indxToPlot(RiskCorrectTrialsIndx),...
                Choice(indxToPlot(RiskCorrectTrialsIndx)),'MarkerFaceColor',[0 0.75 0.25],'MarkerEdgeColor', [0 0.75 0.25]);
            %Plot Risk Lick, punished
            RiskPunishedTrialsIndx   = (OutcomeRecord(indxToPlot) == 3);
            scatter(AxesHandle,  indxToPlot(RiskPunishedTrialsIndx),...
                Choice(indxToPlot(RiskPunishedTrialsIndx)),'MarkerFaceColor','r','MarkerEdgeColor', 'r');
            %Plot Lick, omission
            DidNotChooseTrialsIndx   = (OutcomeRecord(indxToPlot) == 4);
            scatter(AxesHandle,  indxToPlot(DidNotChooseTrialsIndx),...
                Choice(indxToPlot(DidNotChooseTrialsIndx)),'MarkerEdgeColor','b','MarkerFaceColor','b');
            %Plot NoLick, punished
            RandomPunishedTrialsIndx =(OutcomeRecord(indxToPlot)  == 5);
            scatter(AxesHandle,  indxToPlot(RandomPunishedTrialsIndx),...
                Choice(RandomPunishedTrialsIndx),'bo','MarkerFaceColor',[1 1 1]);
            %Plot NoLick, omission
            DidNotChooseTrialsIndx   = (OutcomeRecord(indxToPlot) == 6);
            scatter(AxesHandle,  indxToPlot(DidNotChooseTrialsIndx),...
                Choice(indxToPlot(DidNotChooseTrialsIndx)),'bo','MarkerFaceColor',[1 1 1]);
            xlabel(AxesHandle, 'Trial#', 'FontSize', 18);
            drawnow;
        end

end

end

function [mn,mx] = rescaleX(AxesHandle,CurrentTrial,nTrialsToShow)
 % After this fraction of visible trials, the trial position in the window "sticks" and the window begins to slide through trials.
mn = 1;
mx = nTrialsToShow;
if CurrentTrial > 20
    mn = CurrentTrial - (nTrialsToShow - 1);
    mx = CurrentTrial;
    set(AxesHandle,'XLim',[mn-1 mx+1]);
end
end


