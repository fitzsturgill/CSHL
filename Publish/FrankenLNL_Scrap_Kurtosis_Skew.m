

[cuedData, xData] = phAlignedWindow(TE, trialsByType{1} & hitTrials & ismember(TE.sessionIndex, sessionIndices), 1,...
    'zeroField', 'Us', 'FluorDataField', 'raw', 'PhotometryField', 'Photometry', 'zeroTimes', TE.Us, 'window', window);


data_std = std(cuedData, 0, 1);
data_mean = mean(cuedData);

data_cv = data_std ./ data_mean;
data_fano = data_std.^2;
data_std = smooth(data_std, 5);
data_skew = skewness(cuedData);
data_skew = smooth(data_skew, 5);
data_kurt = kurtosis(cuedData);
data_kurt = smooth(data_kurt, 5);

ensureFigure('test', 1); 
subplot(2,2,1); plot(xData, data_mean);
subplot(2,2,2); plot(xData, data_std.^2);
subplot(2,2,3); plot(xData, data_skew);
subplot(2,2,4); plot(xData, data_kurt);

%%
ensureFigure('test2', 1);
axes;
yyaxis left;
plot(xData, data_mean);
yyaxis right;
plot(xData, data_cv);
grid on;
%% 
ensureFigure('moments');
subplot(2,2,1); plot(xData, data_mean); ylabel('mean');
subplot(2,2,2); plot(xData, data_std.^2); ylabel('variance');
subplot(2,2,3); plot(xData, data_skew); ylabel('skew');
subplot(2,2,4); plot(xData, data_kurt); ylabel('kurtosis');

%% strips
data_mean2 = (data_mean - min(data_mean)) ./ range(data_mean);
cData =[data_mean2; data_std.^2 ./ max(data_std.^2); data_skew ./ max(data_skew)];
ensureFigure('strips_cued', 1); imagesc(cData, [-1 1]);
%%
[uncuedData, xData] = phAlignedWindow(TE, uncuedReward & hitTrials & ismember(TE.sessionIndex, sessionIndices), 1,...
    'zeroField', 'Us', 'FluorDataField', 'raw', 'PhotometryField', 'Photometry', 'zeroTimes', TE.Us, 'window', window);


data_std = std(uncuedData, 0, 1);
data_std = smooth(data_std, 5);
data_mean = mean(uncuedData)';

data_cv = data_std ./ data_mean;
data_fano = data_std.^2;
data_skew = skewness(uncuedData);
data_skew = smooth(data_skew, 5);
data_kurt = kurtosis(uncuedData);

ensureFigure('test', 1); 
subplot(2,2,1); plot(xData, data_cv);
subplot(2,2,2); plot(xData, data_std);
subplot(2,2,3); plot(xData, data_mean);
subplot(2,2,4); plot(xData, data_fano);

%%
ensureFigure('test2', 1);
axes;
yyaxis left;
plot(xData, data_mean);
yyaxis right;
plot(xData, data_fano);
grid on;

%% 
ensureFigure('moments');
subplot(2,2,1); plot(xData, data_mean); ylabel('mean');
subplot(2,2,2); plot(xData, data_std.^2); ylabel('variance');
subplot(2,2,3); plot(xData, data_skew); ylabel('skew');
subplot(2,2,4); plot(xData, data_kurt); ylabel('kurtosis');

%% strips
data_mean2 = (data_mean - min(data_mean)) ./ range(data_mean);
cData =[data_mean2'; data_std'.^2 ./ max(data_std.^2); data_skew' ./ max(data_skew)];
ensureFigure('strips_uncued', 1); imagesc(cData, [-1 1]);