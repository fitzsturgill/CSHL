function port = eventPortFromEventStrings(eventStrings)
    % extracts port number from neuralynx header cell array (see Nlx2MatEV)
    % port -> double vector containing port numbers or NaN for header's
    % without port (e.g. starting recording)
    expr = '(?<=port.).(?=.value)';
    
    port = NaN(size(eventStrings))';
    
    for counter = 1:length(eventStrings)
        es = eventStrings{counter};
        es = regexp(es, expr, 'match'); % should always be scalar
        if ~isempty(es)
            port(counter) = str2double(es);
        end
    end