function settings = bpParseSettings(defaults, inputs)
% point of this function is to combine parameters set by arguments to a
% function call and parameters set by default within the called function.
% settings- ouput is a structure with fields corresponding to the
% initialized parameters (initialized via defaults argument)
% together. That way you can use the output of this function both to supply
% the necessary parameter values iwthin it as well as to STORE the
% parameters used for future reference. In other words, so that you can
% know what you did....
%   ex. defaults = {'param1', 'default1', '';...
%                   'param2', 'default2', ''}
%       inputs   = {'param2', 'usersupplied'};
%    def:   function out = myfunction(defaults, varargin)
%                out = bpParseSettings(defaults, varargin); 
%    function call: s = myfunction(defaults, inputs)    

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
    