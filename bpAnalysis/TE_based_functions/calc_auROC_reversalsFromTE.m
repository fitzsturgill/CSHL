
function auROC = calc_auROC_reversalsFromTE(TE, data, trials1, trials2, window, reset)
% calculate a moving auROC curve comparing data values between two 
% interleaved trial sets (trials1 and trials2) that is reset each reversal
% Assumes that TE contains fields TE.BlockChange and TE.SessionChangee
% parameters: 
% TE- a TE containing BlockChange and SessionChange Fields
% data-  length nTrials data vector
% trials1, trials2, 2 interleaved sets of trials
% window (defualt = 20) moving window for auROC determination
% reset (dault = 1) - whether to reset window at junctions between sessions, blocks

if nargin < 6
    reset = 1;
end

if nargin < 5
    window = 20;
end
data1 = [];
data2 = [];
nTrials = length(TE.filename);
auROC = NaN(nTrials,1);
data = data(:);
for counter = 1:nTrials
    % reset windowed data arrays upon session or block change
    if reset && (counter == 1 || TE.BlockChange(counter) || TE.sessionChange(counter))
        data1 = [];
        data2 = [];
    end
    if length(data1) > window
        data1 = data1(2:end);
    end
    if length(data2) > window
        data2 = data2(2:end);
    end        
    if trials1(counter)
        data1(end+1) = data(counter);
    elseif trials2(counter)
        data2(end+1) = data(counter);
    end
    
    % calc auROC
    auROC(counter) = rocarea(stripNaNs(data1), stripNaNs(data2), 'scale');
end
    