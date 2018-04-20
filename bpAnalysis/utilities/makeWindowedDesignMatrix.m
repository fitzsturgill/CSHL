function [X, truncPoints] = makeWindowedDesignMatrix(x, varargin)
% input: x, unshifted regressor, nSamples x nTrials
% outputs: X, shifted regressor matrix, truncPoints = [startIndex
% stopIndex], start and stop indices to access matching truncated version
% of target/dependent variable

% designed for equal trial sizes...
    %% optional parameters, first set defaults
    defaults = {...
        'window', [-1 1];...
        'Fs', 20;... %
        };
    [s, ~] = parse_args(defaults, varargin{:}); % combine default and passed (via varargin) parameter settings


    % generate weight time offset indices
    Kw2 = s.window * s.Fs; % scaled kernel window
    Kwix = -Kw2(1):-1:-Kw2(2); % offset indices
    
    margins = [sum(Kwix > 0) sum(Kwix < 0)]; % to avoid edge effects due to using circ-shifted versions of the regressor
    kSize = length(Kwix);
    
    nSamples = size(x, 1);
    nTrials = size(x, 2);
    nSamplesTrunc = nSamples - margins(1) - margins(2);
    X = zeros(nSamplesTrunc, nTrials, length(Kwix));
    
    for counter = 1:length(Kwix)
        shifted = circshift(x, Kwix(counter), 1);
        X(:,:,counter) = shifted(margins(1) + 1 : end - margins(2), :);
    end
    
    X = reshape(X, nSamplesTrunc * nTrials, length(Kwix));
    truncPoints = [margins(1) + 1, nSamples - margins(2)];
    
    
    
    