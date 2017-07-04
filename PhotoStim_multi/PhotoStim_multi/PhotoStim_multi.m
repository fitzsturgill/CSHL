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
function PhotoStim_multi
% PhotoStim w/ multiple frequencies for optogenetic tagging
% Written by Hyun-Jae Pi 11/2014.

global BpodSystem
PulsePal;

%% Define parameters
S = BpodSystem.ProtocolSettings; % Load settings chosen in launch manager into current workspace as a struct called S

MaxTrials = 160+1;
TrialTypes=repmat([1 2 3 4 5],1,60);
BpodSystem.Data.TrialTypes = []; % The trial type of each trial completed will be added here.

%% Initialize plots
BpodSystem.ProtocolFigures.OutcomePlotFig = figure('Position', [400 600 1000 200],'Name','Outcome plot','numbertitle','off', 'MenuBar', 'none', 'Resize', 'off');
BpodSystem.GUIHandles.OutcomePlot = axes('Position', [.075 .3 .89 .6]);
OutcomePlot_PhotoStim_multi(BpodSystem.GUIHandles.OutcomePlot,'init',2-TrialTypes,MaxTrials-1);

%% Main trial loop
for currentTrial = 1:MaxTrials-1
    switch TrialTypes(currentTrial) % Determine trial-specific state matrix fields
        case 1
           load('D:\Data\BPOD_Quentin\Protocols\PhotoStim_multi\LightTrain_5hz_1ms.mat');
        case 2
           load('D:\Data\BPOD_Quentin\Protocols\PhotoStim_multi\LightTrain_10hz_1ms.mat');
        case 3
           load('D:\Data\BPOD_Quentin\Protocols\PhotoStim_multi\LightTrain_20hz_1ms.mat');
        case 4
           load('D:\Data\BPOD_Quentin\Protocols\PhotoStim_multi\LightTrain_40hz_1ms.mat');
        case 5
           load('D:\Data\BPOD_Quentin\Protocols\PhotoStim_multi\LightTrain_80hz_1ms.mat');
    end
    ProgramPulsePal(ParameterMatrix);
    
   % S = BpodParameterGUI('sync', S); % Sync parameters with BpodParameterGUI plugin
    
    sma = NewStateMatrix(); % Assemble state matrix  
    if TrialTypes(currentTrial)==1
    sma = AddState(sma, 'Name', 'DeliverStimulus', ...
        'Timer', 0,...
        'StateChangeConditions', {'Tup', 'LightTrain_5hz'},...
        'OutputActions', {});
    elseif TrialTypes(currentTrial)==2
    sma = AddState(sma, 'Name', 'DeliverStimulus', ...
        'Timer', 0,...
        'StateChangeConditions', {'Tup', 'LightTrain_10hz'},...
        'OutputActions', {});
    elseif TrialTypes(currentTrial)==3
    sma = AddState(sma, 'Name', 'DeliverStimulus', ...
        'Timer', 0,...
        'StateChangeConditions', {'Tup', 'LightTrain_20hz'},...
        'OutputActions', {});   
    elseif TrialTypes(currentTrial)==4
    sma = AddState(sma, 'Name', 'DeliverStimulus', ...
        'Timer', 0,...
        'StateChangeConditions', {'Tup', 'LightTrain_40hz'},...
        'OutputActions', {});     
    else
    sma = AddState(sma, 'Name', 'DeliverStimulus', ... 
        'Timer', 0,...
        'StateChangeConditions', {'Tup', 'LightTrain_80hz'},...
        'OutputActions', {});  
    end
    sma = AddState(sma, 'Name', 'LightTrain_5hz', ...
        'Timer', 0,...
        'StateChangeConditions', {'Tup', 'ITI'},...
        'OutputActions', {'BNCState',2});
    sma = AddState(sma, 'Name', 'LightTrain_10hz', ...
        'Timer', 0,...
        'StateChangeConditions', {'Tup', 'ITI'},...
        'OutputActions', {'BNCState',2});
    sma = AddState(sma, 'Name', 'LightTrain_20hz', ...
        'Timer', 0,...
        'StateChangeConditions', {'Tup', 'ITI'},...
        'OutputActions', {'BNCState',2});
    sma = AddState(sma, 'Name', 'LightTrain_40hz', ...
        'Timer', 0,...
        'StateChangeConditions', {'Tup', 'ITI'},...
        'OutputActions', {'BNCState',2});
    sma = AddState(sma, 'Name', 'LightTrain_80hz', ...
        'Timer', 0,...
        'StateChangeConditions', {'Tup', 'ITI'},...
        'OutputActions', {'BNCState',2});
    sma = AddState(sma, 'Name', 'ITI', ...
        'Timer', 10,...
        'StateChangeConditions', {'Tup', 'exit'},...
        'OutputActions', {});  
    
    SendStateMatrix(sma);
    RawEvents = RunStateMatrix;
    if ~isempty(fieldnames(RawEvents)) % If trial data was returned
        BpodSystem.Data = AddTrialEvents(BpodSystem.Data,RawEvents); % Computes trial events from raw data
        % BpodSystem.Data = BpodNotebook(BpodSystem.Data); % Sync with Bpod notebook plugin
        BpodSystem.Data.TrialSettings(currentTrial) = S; % Adds the settings used for the current trial to the Data struct (to be saved after the trial ends)
        BpodSystem.Data.TrialTypes(currentTrial) = TrialTypes(currentTrial); % Adds the trial type of the current trial to data
        UpdateOutcomePlot(TrialTypes, BpodSystem.Data);
        SaveBpodSessionData; % Saves the field BpodSystem.Data to the current data file
    end
    
%     if BpodSystem.Status.BeingUsed == 0
%         return
%     end
    
    % close the procotol
    if currentTrial == MaxTrials-1
        BpodSystem.Status.BeingUsed = 1;
        RunProtocol;
    end
    
end

function UpdateOutcomePlot(TrialTypes, Data)
global BpodSystem
Outcomes = zeros(1,Data.nTrials);
for x = 1:Data.nTrials
    if ~isnan(Data.RawEvents.Trial{x}.States.DeliverStimulus(1))
        Outcomes(x) = 1;
    else
        Outcomes(x) = 3;
    end
end
OutcomePlot_PhotoStim_multi(BpodSystem.GUIHandles.OutcomePlot,'update',Data.nTrials+1,2-TrialTypes,Outcomes)
