function [out1 out2] = calc_selectivity_OutcAuditoryPM(cellid,EpochName,nboot,eventtype)
% 
% CALC_SELECTIVITY_SIDE
%
% Calculates the area under the ROC measure for side selectivity in the
% requested epoch and a bootsrapped P-value.
%
% calc_selectivity_side(cellid,EpochName,nboot)

% AK 10/06 SPR 2010-08-04 PM 02/17

if nargin<4,
    eventtype='behav';
end
switch eventtype
    case 'stim'
        ST = loadcb(cellid,'STIMSPIKES');
        TE = loadcb(cellid,'StimEvents');
    case 'behav'
        ST = loadcb(cellid,'EVENTSPIKES');
        TE = loadcb(cellid,'TrialEvents');
end

if strcmpi(cellid,'default')
    out = [strcat({'Dside_','Pside_'},{EpochName,EpochName})];
    return;
end

epoch_pos = findcellstr(ST.epochs(:,1),EpochName);

if (epoch_pos == 0)
    error('Epoch name not found');
end

RATE = ST.epoch_rates{epoch_pos};

% try
    ValidTrials = selecttrial(TE,'CompletedTrial==1 & Modality==2');
    CorrectAudTrials        = intersect(find(TE.CorrectChoice==1),ValidTrials);
    WrongAudTrials       = intersect(find(TE.CorrectChoice==0),ValidTrials);
% catch
%     ValidTrials=1:length(TE.TrialStart);
%     posLEFT = intersect(find(TE.SoundID==1),ValidTrials);
%     posRIGHT = intersect(find(TE.SoundID==2),ValidTrials);
% end

[Doutmod, Poutmod] = rocarea(RATE(CorrectAudTrials),RATE(WrongAudTrials),'boot',nboot,'scale');

out1 =  Doutmod ;
out2 = Poutmod;
