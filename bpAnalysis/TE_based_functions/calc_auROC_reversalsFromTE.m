
function varargout = calc_auROC_reversalsFromTE(TE, data, trials1, trials2, varargin)% window, reset)
% calculate a moving auROC curve comparing data values between two 
% interleaved trial sets (trials1 and trials2) that is reset each reversal
% Assumes that TE contains fields TE.BlockChange and TE.SessionChangee
% INPUT PARAMETERS: 
% TE- a TE containing BlockChange and SessionChange Fields
% data-  length nTrials data vector
% trials1, trials2, 2 interleaved sets of trials
% window (defualt = 20) moving window for auROC determination
% reset (dault = 1) - whether to reset window at junctions between sessions, blocks
% OUTPUTS: 
% 1) auROC value, scaled between -1 to 1
% 2) p value for auROC significance if nBoot > 0 
% 3) 95% confidence intervals for auROC value, scaled between -1 to 1

defaults = {...
    'window', 20;...
    'windowMode', 'local';...  % 'local' or 'global', local: distributions gathered from last N instances of "trials1" and "trials2" global: distributions from within last N total trials
    'reset', 1;...
    'nBoot', 0;...
    };
[s, ~] = parse_args(defaults, varargin{:}); % combine default and passed (via varargin) parameter settings

data1 = [];
data2 = [];
nTrials = length(TE.filename);
auROC = NaN(nTrials,1);
pROC = NaN(nTrials, 1);
ciROC = NaN(nTrials, 2);
data = data(:);
for counter = 1:nTrials
    % reset windowed data arrays upon session or block change
    if s.reset && (counter == 1 || TE.BlockChange(counter) || TE.sessionChange(counter))
        data1 = []; % local mode
        data2 = []; % local mode
        lastReset = counter; % global mode
    end
    switch s.windowMode
        case 'local'
            if length(data1) > s.window
                data1 = data1(2:end);
            end
            if length(data2) > s.window
                data2 = data2(2:end);
            end        
            if trials1(counter)
                data1(end+1) = data(counter);
            elseif trials2(counter)
                data2(end+1) = data(counter);
            end
        case 'global'
            startWindow = max(max(1, counter - s.window + 1), lastReset);
            thisWindow = data(startWindow:counter);
            data1 = thisWindow(trials1(startWindow:counter));
            data2 = thisWindow(trials2(startWindow:counter));            
    end
    %     calc auROC, normal version (with hypothesis test)
    if s.nBoot
        [D, P, CI] = rocarea_CI(stripNaNs(data1), stripNaNs(data2), 'boot', s.nBoot, 'scale');
        pROC(counter) = P;
        ciROC(counter, :) = CI;
    else
        D = rocarea_CI(stripNaNs(data1), stripNaNs(data2), 'scale');
    end
    auROC(counter) = D;
%     %% debugging kludge
%     if ismember(counter, [494 497 498])
%         varargout{1} = data1;
%         varargout{2} = data2;
%         varargout{3} = [];
%         return
%     end
end

if nargout > 0
    varargout{1} = auROC;
end
if nargout > 1
    varargout{2} = pROC;
end
if nargout > 2
    varargout{3} = ciROC;
end
    