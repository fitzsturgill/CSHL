function port = eventPortFromEventStrings(eventStrings)

    expr = '(?<=port.).(?=.value)';
    
    port = zeros(size(eventStrings));
    
    for counter = 1:length(eventStrings)
        es = eventStrings{counter};
        portNumber = str2num(regexp(es, expr, 'match')); % should always be scalar
        port(counter) = portNumber;
    end