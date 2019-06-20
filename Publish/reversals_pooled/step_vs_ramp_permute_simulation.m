
% make random steps


nRev = 100;
nTrials = 100;
revPoint = 30;
spread = 20;
stepTrials = round(randn(nRev * 2, 1) * spread/2 + revPoint + spread);
stepTrials = stepTrials(stepTrials > revPoint);
stepTrials = stepTrials(1:nRev);

data = zeros(nRev, nTrials);

for counter = 1:nRev
    data(counter, stepTrials(counter):end) = 1;
end

% permute the data across reversals

% generate subscript indices for permutation
row_ix = repmat(1:nTrials, nRev, 1);
col_ix = zeros(nRev, nTrials);
for counter = 1:nTrials
    col_ix(:,counter) = randperm(sum(nRev))';
end

lin_ix = sub2ind([nRev nTrials], col_ix, row_ix);
perm_data = data(lin_ix);
perm_data_smooth = smoothdata(perm_data, 2, 'movmean');

% do changepoint detection
tic
data_cp = bpChangePoints(data, 2, 10, 'up');
toc
disp('first');
tic
perm_data_cp = bpChangePoints(perm_data, 2, 10, 'up');
disp('second');
toc
%
[data_aligned, xData] = alignedDataWindow(data, 1:nRev, 'zeroTimes', data_cp.index, 'window', [-spread * 2 spread * 2], 'Fs', 1, 'startTimes', repmat(1, nRev, 1));
[perm_data_aligned, xData] = alignedDataWindow(perm_data, 1:nRev, 'zeroTimes', perm_data_cp.index, 'window', [-spread * 2 spread * 2], 'Fs', 1, 'startTimes', repmat(1, nRev, 1));


ensureFigure('step_vs_ramp_permute_simulation', 1); 

subplot(3,2,1); imagesc(data); title('non-permuted');
subplot(3,2,3); imagesc(data_aligned); 
subplot(3,2,2); imagesc(perm_data); title('permuted');
subplot(3,2,4); imagesc(perm_data_aligned);


subplot(3,2,5); plot(nanmean(data_aligned)); hold on; 
plot(nanmean(perm_data_aligned)); legend({'non-permuted', 'permuted'}, 'Location', 'best'); legend boxoff;
subplot(3,2,6); colorbar;

formatFigurePoster([10 6], gcf, 12);
