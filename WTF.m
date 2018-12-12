%% fit exponentials to post-reversal data for both acquisition and extinction
% ss = struct('object', zeros(sum(goodReversals), 1), 'gof', zeros(sum(goodReversals), 1), 'output', zeros(sum(goodReversals), 1), 'toFit', zeros(sum(goodReversals), 1));

compFields = {'csPlus', 'csMinus'};
fitFields = {'licks_cs', 'phPeakMean_cs_ch1', 'phPeakMean_cs_ch2'};
outputFields = {'object', 'gof', 'output', 'toFit', 'a', 'b', 'c'};
expFit = struct();
for counter = 1:length(compFields)
    for counter2 = 1:length(fitFields)
        for counter3 = 1:length(outputFields)
            if ~ismember(outputFields{counter3}, {'a', 'b', 'c'})
                expFit.(compFields{counter}).(fitFields{counter2}).(outputFields{counter3}) = cell(sum(goodReversals), 1);
            else
                expFit.(compFields{counter}).(fitFields{counter2}).(outputFields{counter3}) = NaN(sum(goodReversals), 1);
            end
        end
    end
end

expModel = 'a + b * exp(-x/c)';
% these options should work for newCsPlus and newCsMinus (up or down and
% for z scored and lick rate data (both should fall within interval of -100
% to 100
fo = fitoptions('Method', 'NonlinearLeastSquares', 'Robust', 'On',...
    'Upper', [100  100 1000],...
    'Lower', [-100 -100 0],...    
    'StartPoint', [0 1 10]...
    );
for compCounter = 1:length(compFields)
    for fieldCounter = 1:length(fitFields)
        expFitData = AR.(compFields{compCounter}).(fitFields{fieldCounter}).after(goodReversals, :);
        for counter = 1:size(expFitData, 1)        
            toFit = expFitData(counter, ~isnan(expFitData(counter, :)));
            ft = fittype(expModel, 'options', fo);
            try
                [fitobject, gof, output] = fit((0:length(toFit) - 1)', toFit', ft, fo);
                expFit.(compFields{compCounter}).(fitFields{fieldCounter}).object{counter} = fitobject;
                expFit.(compFields{compCounter}).(fitFields{fieldCounter}).a(counter) = fitobject.a;
                expFit.(compFields{compCounter}).(fitFields{fieldCounter}).b(counter) = fitobject.b;
                expFit.(compFields{compCounter}).(fitFields{fieldCounter}).c(counter) = fitobject.c;
            catch
                continue
            end

        end
    end
end
