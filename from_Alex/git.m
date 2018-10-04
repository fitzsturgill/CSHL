function git(varargin)
%
% Passes on all inputs to the shell git function.
%
% Inputs  :: varargin of strings that will be passed on to shell git via system()
%            Last input may be '--test' or '--dry-run'.  If either is the case, 
%            this function runs in test  mode (ie., --dry-run mode for git).
%
% Outputs :: None.  Will display the return value from git.
%
% Comments : Yes, this ~50 lines of code allows you to avoid just typing "!git" instead.  You're welcome.
%
% Alex Vaughan, 2015

verbose = 0;

% Parse last input for --dry-run variable.
test_only = 0;
if strcmp(varargin{end},'--test') || strcmp(varargin{end},'--dry-run')
    test_only = 1;
    inputs = {varargin{1:end-1}};
else
    inputs = varargin;
end

% Eval all inputs and make sure strings are still escaped.
git_strings = cell(size(inputs));
for i = 1:length(inputs)
    if verbose,
        fprintf('inputs{%.0f} :: %s\n',i,inputs{i})
    end
    if isempty(strfind(inputs{i},' '))
        git_strings{i} = inputs{i};
    else
        git_strings{i} =  sprintf('''%s''',inputs{i});
    end
end

% Assemble and run git command.
git_command = strjoin({'git',git_strings{:}});
if test_only
    fprintf('git.m :: Running TEST only.\n')
    git_command = strjoin({git_command,'--dry-run'});
end

[~,cwd,~] = fileparts(pwd);
fprintf('git.m :: In folder %s :: %s\n',cwd,pwd);
fprintf('git.m :: Running command ::\n\t.\t%s\n',git_command)
[status,cmdout] = system(git_command);


% Print output of command.
cmdout = strsplit(cmdout,'\n');
cmdout = strjoin({'',cmdout{:}},'\n\t.\t');
fprintf('git.m :: status :: %.0f\n',status)
fprintf('git.m :: cmdout :: %s\n',cmdout)
%fprintf('git.m :: system() cmdout :: \n"""""<git cmdout>"""""%s\n""""</git cmdout>"""""\n',cmdout)
