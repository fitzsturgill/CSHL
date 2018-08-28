function parsed_params = parse_defaults(default_params,override_params)
% Evaluates two structs and returns a struct similar to "default"
% with specific values altered by the contents of override
%
% Inputs :
%   default_params :: cell or struct.
%   override_paramas :: cell or struct.
%
%   If either input is a cell, must be of format
%   {'field1',value1,'field2',value2, ...}
%
% Outputs :
%   parsed_params :: a struct with the same fields as default_params,
%                    with values altered to match any values found in
%                   override_params
%
%                   This is typically processed in the caller function as
%                       fn = fieldnames(parsed_params)
%                       for i = 1:length(fn)
%                           eval(sprintf('%s = parsed_params.%s;',fn{i},fn{i});
%                       end
%
% Alex Vaughan, 2015

% Convert first input into cell format.
switch class(default_params)
    case 'struct'
        % Convert to cell format, keeping field names
        default_params = struct2cell_withfields(default_params);
    case 'cell'
        % do nothing
    otherwise
        error('Invalid input - must be struct or cell')
end

% Convert second input into cell format.
switch class(override_params)
    case 'struct'
        % Convert to cell format, keeping field names
        override_params = struct2cell_withfields(override_params);
    case 'cell'
        % do nothing
    otherwise
        error('Invalid input - must be struct or cell')
end

if length(override_params)
    
    p = inputParser();
    for d = 1:2:length(default_params)
        addParameter(p,default_params{d},default_params{d+1});
    end
    
    % Do the actual parsing.
    parse(p,override_params{:})
    
    % Return values as parsed_params.whatever
    for d = 1:2:length(default_params)
        eval(sprintf('parsed_params.%s = p.Results.%s;',default_params{d},default_params{d}));
    end
    
else   
    
    % Just turn default_params into a struct;
    for d = 1:2:length(default_params)
        parsed_params.(default_params{d}) = default_params{d+1};
    end
    
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function output_cell = struct2cell_withfields(input_struct)
% Transforms a struct into a cell with format
% {'field1',value1,'field2',value2, ...}
%
% Input : a struct
% Output: a cell, suitably transformed.
%
% Alex Vaughan, 2015

assert(isstruct(input_struct))

fields = fieldnames(input_struct);
values = struct2cell(input_struct);
output_cell = cell(length(fields)*2,1);
for f = 1:length(fields)
    output_cell{(2*(f-1))+1} = fields{f};
    output_cell{(2*(f-1))+2} = values{f};
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%