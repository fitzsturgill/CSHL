%{ 
dev, NMF deconvolution
FINAL SETTINGS: use Fs = 20,  deconvType = 'OASIS', sensorTau = 0.2
%}
importRepositories({'Suite2P', 'OASIS_matlab'});
load('Z:\SummaryAnalyses\wheel_v1\ACh_3_wheel_v1_Feb04_2019_Session1\TE.mat');

%%

photometryField = 'Photometry';


%%
data = TE.(photometryField).data(1).ZS';
data2 = TE.(photometryField).data(2).ZS';


F = [data(:) data2(:)];

%%
ops.fs = TE.(photometryField).sampleRate;
ops.sensorTau = 0.2;
ops.estimateNeuropil = 0;
ops.deconvType = 'OASIS'; % L0 or OASIS


% looks like I need to properly install suite2p...
[sp, ca, coefs, B, sd, ops, baselines] = wrapperDECONV(ops, F);


%% graphs

ensureFigure('test_deconv_scatter', 1); scatter(sp(:,1), sp(:,2))

ensureFigure('test_deconv_hist', 1); histogram(sp(:,1), 0:1:100); hold on; histogram(sp(:,2), 0:1:100); set(gca, 'YScale', 'log');

%% add reconstructed data to TE
sz = size(TE.(photometryField).data(1).ZS);
TE.(photometryField).data(1).reconstructed = reshape(ca(:,1), [sz(2) sz(1)])';
TE.(photometryField).data(2).reconstructed = reshape(ca(:,2), [sz(2) sz(1)])';
TE.(photometryField).data(1).spikes = reshape(sp(:,1), [sz(2) sz(1)])';
TE.(photometryField).data(2).spikes = reshape(sp(:,2), [sz(2) sz(1)])';

%%

trial = 3;
xlim = [5 20];
ensureFigure('test', 1); 
subplot(2,1,1); plot(TE.(photometryField).xData, [zscore(TE.(photometryField).data(1).ZS(trial,:))' zscore(TE.(photometryField).data(1).reconstructed(trial,:))'])
yyaxis right; stem(TE.(photometryField).xData, TE.(photometryField).data(1).spikes(trial,:));
set(gca, 'XLim', xlim);
subplot(2,1,2); plot(TE.(photometryField).xData, [zscore(TE.(photometryField).data(2).ZS(trial,:))' zscore(TE.(photometryField).data(2).reconstructed(trial,:))'])
yyaxis right; stem(TE.(photometryField).xData, TE.(photometryField).data(2).spikes(trial,:));
set(gca, 'XLim', xlim);


%% find optimal tau

taus = 0.01:0.01:2;
ops.fs = TE.(photometryField).sampleRate;
ops.estimateNeuropil = 0;
ops.deconvType = 'OASIS';

data = TE.(photometryField).data(1).ZS';
data2 = TE.(photometryField).data(2).ZS';

F = [data(:) data2(:)];

corrValues = zeros(length(taus), 2);
corrValues2 = zeros(length(taus), 2);
for counter = 1:length(taus)    
    disp(num2str(counter/length(taus)));
    ops.sensorTau = taus(counter);
    % looks like I need to properly install suite2p...
    [sp, ca, coefs, B, sd, ops, baselines] = wrapperDECONV(ops, F);
    
    corrValues(counter, 1) = corr(zscore(ca(:,1)), zscore(data(:)));
    corrValues(counter, 2) = corr(zscore(ca(:,2)), zscore(data(:)));
    
    corrValues2(counter, 1) = corr(zscore(ca(:,1)), zscore(data(:)), 'Type', 'Spearman');
    corrValues2(counter, 2) = corr(zscore(ca(:,2)), zscore(data(:)), 'Type', 'Spearman');    
end

%% plot correlation vs sensor tau

ensureFigure('Pearson_vs_tau', 1); plot(taus, corrValues); xlabel('sensor tau (s)'); title('Pearson'); ylabel('rho'); legend({'ch1', 'ch2'}); legend boxoff;

ensureFigure('Spearman_vs_tau', 1); plot(taus, corrValues2); xlabel('sensor tau (s)'); title('Spearman'); ylabel('rho'); legend({'ch1', 'ch2'}); legend boxoff;