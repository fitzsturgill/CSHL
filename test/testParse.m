function g = testParse(varargin)

% add default settings
    defaults = {...
        'nChannels', 2, '';... % parameter, default value, validation function
        'blStart', 1, ''...
        };
%%
%     prs = inputParser;
%     addParamValue(prs,'nChannels',2)   % output argument names
%     addParamValue(prs,'blStart',1)   % mandatory input arguments
%     parse(prs, varargin{:})
%     g = prs.Results;
  %%
    g = bpParseSettings(defaults, varargin);
end

function settings = bpParseSettings(defaults, inputs)
    
%   ex. defaults = {'param1', 'default1', '';...
%                   'param2', 'default2', ''}
%       inputs   = {'param2', 'usersupplied'};
%       function out = myfunction(defaults, varargin)
%                out = 
%        
% s = bpParseSettings(defaults, varargin); 
% in this example- varargin is the cell array in enclosing function
% containing parameter value pairs, see varargin documentation
% defaults is a n x 3 cell array- each row defines defaults for a given parameter
% First column is a parameter name,
% Second column is a default value
% Third column elements are either empty ('' or [], no validation) or contain a
% validation function handle

    prs = inputParser;
    for i = 1:size(defaults, 1)
        if ~isempty(defaults{i, 3})
            addParameter(prs, defaults{i, 1}, defaults{i, 2}, defaults{i, 3}); % 3rd element validation function handle
        else
            addParameter(prs, defaults{i, 1}, defaults{i, 2}); % no validation        
        end
    end
    parse(prs, inputs{:});
    settings = prs.Results;
end
    