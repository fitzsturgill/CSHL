function list_variables(varargin)
% Takes all variables in the current workspace and lists each variable and
% its contents:

% Indentation for nested variables.
st = dbstack;
stack_depth = -1;
for i = 1:length(st)
    if strcmp(st(i).name,'list_variables')
        stack_depth = stack_depth + 1;
    end
end
%fprintf('\nstack_depth :: %.0f\n',stack_depth)

%  Header and gather variables to list.
indent_prefix = '';
if nargin == 0
    % If nargin == 0, then we list ALL variables in the calling namespace.
    fprintf('\n*** All variables in current workspace. ***\n')
    all_variable_names = evalin('caller','who');
    n_variables = length(all_variable_names);
    all_variable_data = cell(n_variables,1);
    for i = 1:n_variables,
        all_variable_data{i} = evalin('caller',all_variable_names{i});
    end
elseif (stack_depth > 0) && (nargin == 1) && isstruct(varargin{1}),
    % We have been passed a single struct, so unwrap it as if it was a list of
    % variables.
    all_variable_names = fieldnames(varargin{1});
    all_variable_data  = struct2cell(varargin{1});
    n_variables = length(all_variable_names);
else
    % Otherwise we take the inputs and grab the name from the calling space
    % as well as the value from the passed function.4
    if stack_depth == 0,
        fprintf('*** Listing %.0f passed variables. ***\n',nargin)
    end
    n_variables = nargin;
    [all_variable_names,all_variable_data] = deal(cell(n_variables,1));
    for i = 1:n_variables,
        all_variable_names{i} = inputname(i);
        all_variable_data{i} = varargin{i};
    end
        
end

% Set up indentation
max_name_length = 0;
for i = 1:n_variables
    max_name_length = max(max_name_length,length(all_variable_names{i}));
end
max_name_length = max(max_name_length,30);
max_name_length = max_name_length - 2*stack_depth;

% Set up prefix for nested structs
if stack_depth,
    indent_prefix = [repmat([' '],1,(stack_depth)*(3),1),'.'];
else
    % Depth = 0, so don't indent.  Also, put a header.
    indent_prefix = '';
    fprintf(repmat('-',100,1)); fprintf('\n')
    fprintf('Variable Name %s Class \t Size \t\t Value\n',repmat(' ',max_name_length-10,1))
    fprintf(repmat('-',100,1)); fprintf('\n')
end




% Iterate over variables
for i = 1:n_variables
    
    this_variable = all_variable_data{i};
    this_variable_name = all_variable_names{i};    
    var_name_length = length(all_variable_names{i});
    
    % Spaces or row markers.
    spaces_length = max_name_length - var_name_length + 2;
    if mod(i,5)
        spaces_before_value = repmat([' '],1,spaces_length);
    else
        spaces_before_value = repmat(['.'],1,spaces_length);

%        spaces_before_value = repmat(['- '],1,floor(spaces_length/2)+mod(spaces_length,2));
    end
    
    % Print appropriate output for each variable type
    switch class(this_variable)
        case {'double','single','int8','uint8','int16','uint16','int32','uint32','int64','uint64'}
            v_size = size(this_variable);
            if sum(v_size(:)) > 11
                fprintf('%s%s%s:: %s \t size(%s)\n'                ,indent_prefix, this_variable_name, spaces_before_value, class(this_variable), mat2str(size(this_variable)) );
            else
                fprintf('%s%s%s:: %s \t size(%s) \t value: %s\n'   ,indent_prefix, this_variable_name, spaces_before_value, class(this_variable), mat2str(size(this_variable)), mat2str(this_variable,3));
            end
        case 'cell'
            fprintf('%s%s%s:: %s \t size(%s)\n'                 ,indent_prefix, this_variable_name, spaces_before_value, class(this_variable), mat2str(size(this_variable)) );
        case 'char'
            fprintf('%s%s%s:: %s \t''%s''\n'                    ,indent_prefix, this_variable_name, spaces_before_value, class(this_variable), this_variable );
        case 'logical'
            v_size = size(this_variable);
            if sum(v_size(:)) > 2
                fprintf('%s%s%s:: %s \t size(%s)\n'             ,indent_prefix, this_variable_name, spaces_before_value, class(this_variable), mat2str(size(this_variable)) );
            else
                true_false = {'FALSE','TRUE'};
                fprintf('%s%s%s:: %s \t\t\t %s\n'                   ,indent_prefix, this_variable_name, spaces_before_value, class(this_variable), true_false{this_variable+1});
            end
        case 'struct'
            fprintf('%s%s%s:: struct \t size(%s)\n',...
                indent_prefix,...
                this_variable_name,...
                spaces_before_value,...
                mat2str(size(this_variable)));
            if prod(size(this_variable)) == 1,
                list_variables(this_variable);
            end
        otherwise
            fprintf('%s%s %s :: %s\n',indent_prefix,all_variable_names{i},spaces_before_value , class(this_variable));
    end
end