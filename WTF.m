compFields = {'csPlus', 'csMinus';...  % first row is after reversal
              'csMinus', 'csPlus'...   % second row is before reversal
              };
fitFields = {'licks_cs', 'phPeakMean_cs_ch1', 'phPeakMean_cs_ch2'};
outputFields = {'object', 'gof', 'output', 'toFit', 'a', 'b', 'c', 'd', 'tt20'};
weibull = struct();

nReversals = size(AR.csPlus.globalTrialNumber.after, 1);

for counter = 1:length(compFields)
    for counter2 = 1:length(fitFields)
        for counter3 = 1:length(outputFields)
            if ~ismember(outputFields{counter3}, {'a', 'b', 'c', 'd', 'tt20'})
                weibull.(compFields{counter}).(fitFields{counter2}).(outputFields{counter3}) = cell(nReversals, 1);
            else
                weibull.(compFields{counter}).(fitFields{counter2}).(outputFields{counter3}) = NaN(nReversals, 1);
            end
        end
    end
end

weibullModel =  'a * (1 - exp(-1 * (x/b)^c)) + d'; % weibull function, CDF form

maxLickRate = 10;
for compCounter = 1:size(compFields, 2)    
    for fieldCounter = 1:length(fitFields)
        weibullData = [AR.(compFields{2, compCounter}).(fitFields{fieldCounter}).before(:, end - baselineTrials + 1:end) AR.(compFields{1, compCounter}).(fitFields{fieldCounter}).after];
        for counter = 1:size(weibullData, 1)        
            toFit = weibullData(counter, ~isnan(weibullData(counter, :)));
            switch compFields{1, compCounter}
                case 'csPlus'
                    fo = fitoptions('Method', 'NonlinearLeastSquares',... 
                        'Upper', [maxLickRate  Inf Inf maxLickRate],...  % 20 (3rd upper)
                        'Lower', [0 0 0 0],...    % 'Lower', [0 0 -1/5 0 -1/5],...                    
                        'StartPoint', [range(toFit) baselineTrials baselineTrials min(toFit)]...  % 'StartPoint', [range(toFit) baselineTrials baselineTrials min(toFit)]...
                        );
                case 'csMinus'
                    fo = fitoptions('Method', 'NonlinearLeastSquares',... 
                        'Upper', [0  Inf Inf maxLickRate],...  % 20 (3rd upper)
                        'Lower', [-maxLickRate 0 0 0],...    % 'Lower', [0 0 -1/5 0 -1/5],...                    
                        'StartPoint', [-range(toFit) baselineTrials baselineTrials max(toFit)]... % 'StartPoint', [-range(toFit) baselineTrials baselineTrials max(toFit)]...
                        );         
                    
%                 case 'csPlus'
%                     fo = fitoptions('Method', 'NonlinearLeastSquares',... 
%                         'Upper', [maxLickRate  Inf Inf maxLickRate],...  % 20 (3rd upper)
%                         'Lower', [0 0 0 0],...    % 'Lower', [0 0 -1/5 0 -1/5],...                    
%                         'StartPoint', [range(toFit) baselineTrials baselineTrials min(toFit)]...  % 'StartPoint', [range(toFit) baselineTrials baselineTrials min(toFit)]...
%                         );
%                 case 'csMinus'
%                     fo = fitoptions('Method', 'NonlinearLeastSquares',... 
%                         'Upper', [0  Inf Inf maxLickRate],...  % 20 (3rd upper)
%                         'Lower', [-maxLickRate 0 0 0],...    % 'Lower', [0 0 -1/5 0 -1/5],...                    
%                         'StartPoint', [-range(toFit) baselineTrials baselineTrials max(toFit)]... % 'StartPoint', [-range(toFit) baselineTrials baselineTrials max(toFit)]...
%                         );                       
            end            
            ft = fittype(weibullModel, 'options', fo);
            weibull.(compFields{1, compCounter}).(fitFields{fieldCounter}).toFit{counter} = toFit;
            try
                [fitobject, gof, output] = fit((0:length(toFit) - 1)', toFit', ft, fo);
                weibull.(compFields{1, compCounter}).(fitFields{fieldCounter}).object{counter} = fitobject;
                weibull.(compFields{1, compCounter}).(fitFields{fieldCounter}).a(counter) = fitobject.a;
                weibull.(compFields{1, compCounter}).(fitFields{fieldCounter}).b(counter) = fitobject.b;
                weibull.(compFields{1, compCounter}).(fitFields{fieldCounter}).c(counter) = fitobject.c;
                weibull.(compFields{1, compCounter}).(fitFields{fieldCounter}).d(counter) = fitobject.d;
                weibull.(compFields{1, compCounter}).(fitFields{fieldCounter}).output{counter} = output;
                weibull.(compFields{1, compCounter}).(fitFields{fieldCounter}).gof{counter} = gof;  
                if counter == 3
                    disp('wtf');
                end
                switch compFields{1, compCounter}
                    case 'csPlus'
                        thresh = fitobject(0) + (fitobject(Inf) - fitobject(0)) * 0.2;
                        latency = find(toFit >= thresh, 1) - baselineTrials;
                        if isempty(latency)
                            continue
                        end
                        weibull.(compFields{1, compCounter}).(fitFields{fieldCounter}).tt20(counter) = latency;
                    case 'csMinus'
                        thresh = fitobject(0) - (fitobject(Inf) - fitobject(0)) * 0.2;
                        latency = find(toFit <= thresh, 1) - baselineTrials;
                        if isempty(latency)
                            continue
                        end
                        weibull.(compFields{1, compCounter}).(fitFields{fieldCounter}).tt20(counter) = latency;
                end                
            catch
                continue
            end

        end
    end
end

%%

% dT = 0.05;
% 
% sr = 40;
% 
% duration = 60;
% 
% nTrials = 50;
% 
% stimes = cell(nTrials, 1);
% for counter = 1:nTrials
%     stimes{counter} = makeSpikes(dT, sr, duration, 1);
% end
% fitz= 
% binraster = stimes2binraster(stimes, 0:dT:duration, dT, repmat([0 60], nTrials, 1));
% 
% 
% 
% % binraster = stimes2binraster(spikes
% 
% 
% function spikes = makeSpikes(timeStepS, spikesPerS, durationS, numTrains)
% 
% if (nargin < 4)
%     numTrains = 1;
% end
% times = [0:timeStepS:durationS];
% spikes = zeros(numTrains, length(times));
% for train = 1:numTrains
%     vt = rand(size(times));
%     spikes(train, :) = (spikesPerS*timeStepS) > vt;
% end
% end