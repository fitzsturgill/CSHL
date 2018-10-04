function newVarargin = rewrite_varargin(parsed, varargin)
% updates varargin via the output of parse_args in order to pass an updated
% version to functions called by the parent function
% companion function to parse_args

%  parsed, structure outputed by parse_args with contents of structure
%  fields possibly further modified (e.g. to handle certain conditions)

% varargin, just feed in here 'varargin{:}' or 'varargin'

% Fitz Sturgill 7/2018

    % if varargin is passed without expansion
    if length(varargin) == 1 && iscell(varargin{1})
        varargin = varargin{1};
    end
    
    newVarargin = varargin;
    
    fields = fieldnames(parsed);    
    
    for counter = 1:length(fields)
        field = fields{counter};
        wpos = cellfun(@(x) ischar(x) && strcmp(x, field), newVarargin);
        if ~any(wpos)
            newVarargin(end+1:end+2) = {field, parsed.(field)};
        else
            newVarargin{find(wpos)+ 1} = parsed.(field);
        end
    end