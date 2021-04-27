%{
----------------------------------------------------------------------------

This file is part of the Bpod Project
Copyright (C) 2014 Joshua I. Sanders, Cold Spring Harbor Laboratory, NY, USA

----------------------------------------------------------------------------

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, version 3.

This program is distributed  WITHOUT ANY WARRANTY and without even the 
implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
%}
% function OutcomePlot(AxesHandle,TrialTypeSides, OutcomeRecord, CurrentTrial)
function OutcomePlot_Gambling12(AxesHandle, Action, varargin)
%
% My input: OutcomePlot_Gambling12(BpodSystem.GUIHandles.OutcomePlot,'init',1-TrialType, UsOutcome1);
%           OutcomePlot_Gambling12(BpodSystem.GUIHandles.OutcomePlot,'update',Data.nTrials+1,1-TrialType,Outcomes, UsOutcome)%
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
% Based on code provided by Heguedus Panna
%
% Cecília Pardo-Bellver, 2020
% Laboratory of Network Neurophysiology
% Institute of Experimantal Medicine, Hungary.
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
        scatter(AxesHandle,  find(ColorList(1:nTrialsToShow) == 0), TrialTypesSide(ColorList(1:nTrialsToShow) == 0),'MarkerFaceColor',[1 1 1],'MarkerEdgeColor', [0 0 1]);
        hold(AxesHandle, 'on');
        scatter(AxesHandle,  find(ColorList(1:nTrialsToShow) == 1), TrialTypesSide(ColorList(1:nTrialsToShow) == 1),'MarkerFaceColor',[1 1 1],'MarkerEdgeColor', [0 1 0]);
        scatter(AxesHandle,  find(ColorList(1:nTrialsToShow) == 2), TrialTypesSide(ColorList(1:nTrialsToShow) == 2),'MarkerFaceColor',[1 1 1],'MarkerEdgeColor', [1 0 0]);
        set(AxesHandle,'TickDir', 'out','XLim',[0 21],'YLim', [-2, 1], 'YTick', [-1 0],'YTickLabel', {'Type 2 - L','Type 1 - R'}, 'FontSize', 14);
        xlabel(AxesHandle, 'Trial#', 'FontSize', 14);
                
    case 'update'
        CurrentTrial = varargin{1};
        TrialTypesSide = varargin{2};
        OutcomeRecord = varargin{3};
        ColorList = varargin{4};
        
        if CurrentTrial<1
            CurrentTrial = 1;
        end
        
        % recompute xlim
        [mn, mx] = rescaleX(AxesHandle,CurrentTrial,nTrialsToShow);
        
%         %plot future trials
%         FutureTrialsIndx = CurrentTrial;
%         scatter(AxesHandle,  FutureTrialsIndx(ColorList(FutureTrialsIndx) == 0), TrialTypesSide(FutureTrialsIndx(ColorList(FutureTrialsIndx) == 0)),'MarkerFaceColor',[1 1 1],'MarkerEdgeColor', [0 0 1]);
%         hold(AxesHandle, 'on');
%         scatter(AxesHandle,  FutureTrialsIndx(ColorList(FutureTrialsIndx) == 1), TrialTypesSide(FutureTrialsIndx(ColorList(FutureTrialsIndx) == 1)),'MarkerFaceColor',[1 1 1],'MarkerEdgeColor', [0 1 0]);
%         scatter(AxesHandle,  FutureTrialsIndx(ColorList(FutureTrialsIndx) == 2), TrialTypesSide(FutureTrialsIndx(ColorList(FutureTrialsIndx) == 2)),'MarkerFaceColor',[1 1 1],'MarkerEdgeColor', [1 0 0]);
%         
%         %Plot current trial
%         scatter(AxesHandle,CurrentTrial,TrialTypesSide(CurrentTrial), 'o', 'MarkerFaceColor',[1 1 1],'MarkerEdgeColor', 'k')
%         scatter(AxesHandle,CurrentTrial,TrialTypesSide(CurrentTrial), '+', 'MarkerFaceColor',[1 1 1],'MarkerEdgeColor', 'k')
        
        %Plot past trials
        if ~isempty(OutcomeRecord)
            indxToPlot = mn:CurrentTrial-1;
            %Plot NoLick
            EarlyWithdrawalTrialsIndx = (OutcomeRecord(indxToPlot) == 0);
            scatter(AxesHandle,  indxToPlot(EarlyWithdrawalTrialsIndx), TrialTypesSide(indxToPlot(EarlyWithdrawalTrialsIndx)),'bo','MarkerFaceColor','b');
            %Plot Lick L reward
            CorrectTrialsIndx         = (OutcomeRecord(indxToPlot) == 2);
            scatter(AxesHandle,  indxToPlot(CorrectTrialsIndx), TrialTypesSide(indxToPlot(CorrectTrialsIndx)),'MarkerFaceColor','g','MarkerEdgeColor', 'g');
            %Plot Lick R reward
            CorrectTrialsIndx         = (OutcomeRecord(indxToPlot) == 1);
            scatter(AxesHandle,  indxToPlot(CorrectTrialsIndx), TrialTypesSide(indxToPlot(CorrectTrialsIndx)),'MarkerFaceColor','g','MarkerEdgeColor', 'g');
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
