function [g, err] = parse_args(default_args,varargin)
% PARSE_ARGS   Parse input arguments.
%   [G, ERR] = PARSE_ARGS(DEFAULT_ARGS) converts default argument name
%   value pairs to fields of the output structure G. Additional name value
%   pairs can be passed as varargin.
%
%   See also VIEWCELL2B.

%   Edit log: BH 6/23/11, 4/19/12, 7/5/12

% Input argument check; concert to structure
if isstruct(default_args)
    def = default_args;
else
%     d = default_args';  % flip matrix
    d = default_args;  % DON'T flip matrix    
    try
%         def = struct(d{:}); 
        % avoid bug where user supplies a cell array (inadvertant expansion
        % of scalar structure)
        def = convertToStructure(d); % subfunction
    catch ME
        disp(ME.message)
        error('Argument error in default list');
    end
end

% If there are args, convert them into arg list
if ~isempty(varargin)
    try
%         g = struct(varargin{:}); 
        % avoid bug where user supplies a cell array (inadvertant expansion
        % of scalar structure)
        g = convertToStructure(varargin); % subfunction
    catch ME
        disp(ME.message)
        error('Argument error in the {''param'', value} sequence');
    end
else
    g.default_args = 1;
end

% Get all the field names set
g_fields = fieldnames(g);
d_fields = fieldnames(def);

% If fields were set that are not on the default arglist, then these are 
% considered errors
err = setdiff(g_fields,d_fields);

% Fields that weren't set yet
toset_fields = setdiff(d_fields,g_fields);

% Set them to the default values
for iF = 1:length(toset_fields)
    g.(toset_fields{iF}) = def.(toset_fields{iF});
end

function s = convertToStructure(ca) % ca = cell array, e.g. varargin

    % cell array (as supplied by functions that use parseargs) for example 
    % consists of rows of parameter value pairs, we need to make this an alternating vector    
    if size(ca, 1) > 1; 
        ca = ca';
        ca = reshape(ca, 1, numel(ca)); % now it's a row vector
    end
    counter = 1;
    while counter+1 <= length(ca) 
        prop = ca{counter};
        val = ca{counter+1};
        s(1).(prop) = val;  
        counter = counter + 2;
    end