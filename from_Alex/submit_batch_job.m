function [job_outputs,job_ndgrid]= submit_batch_job(function_to_call,function_params,cluster_params)
%
% submit_batch_job() will submit a qsub batch job to the CSHL UGE engine on
% the blackNblue cluster, using the hyperparametersdefined in params_list,
% and the cluster/qsub parameters defined in cluster_params
%
% Inputs ::
%
%   1) function_to_call :: is the name of the function to be passed to qsub.  It is
%       assumed that this is a compiled function accessible to the path.
%
%   2) function_params :: is a cell of cells, carrying of parameters to be passed to
%       function_to_call.  These are interpreted as ranges of
%       hyperparameters with which to call function_to_call, as follows.
%
%       Example usage ::
%           For:
%               function_params = { {1}, {'a','b'} }
%           two jobs will be submitted to execute:
%               function_to_call{1,a)
%               function_to_call{1,b)
%
%   3) cluster_params :: is a struct of default parametes for job
%       submission, including command line flags for qsub.
%
% Outputs ::
%
%   1) job_outputs, a cell array of job outputs.  Each element is simply
%   the encapsulation of all outputs from the job. If function_to_call has multiple
%   outputs, these outputs are encapsulated as a cell.
%
%   2) job_ndgrid, a cell array output from the ndgrid implementation of
%   hyperparameter implementation.  For a job with parameters a,b,c,d, the
%   parameters used to build job_outputs{17} are found via
%   settings are found via job_ndgrid{1}(17), job_ndgrid{2}(17), etc.
%
%
% Alex Vaughan, 2015
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Setup cluster parameters
default_cluster_params = struct(...
 'qsub_command_resources',      '-l num_cpus>15'   ...
);
cluster_params = parse_defaults(default_cluster_params,cluster_params);
cluster_params_fields = fieldnames(default_cluster_params)
for i = 1:length(cluster_params_fields)
    eval('%s = %s',cluster_params_fields{i},default_cluster_params.(default_cluster_params{i}))
end

% Define job ndgrid
n_function_params = length(function_params);
function_params_sizes = zeros(n_function_params,1)
job_ndgrid_inputs = cell(job_ndgrid_inputs,1);
for i = 1:n_function_params
    function_params_sizes(i) = length(function_params{i});
    job_ndgrid_inputs{i} = 1:function_params_sizes(i);
end
[job_ndgrid{:}] = ndgrid(job_ndgrid_inputs{:});
job_ndgrid_length = numel(job_ndgrid{1});


for nd = 1:job_ndgrid_length
    
    for np = n_function_params,
        if isnumeric(function_params{nd})
            this_function_params{np} = function_params{nd}(job_ndgrid{np}(nd));
        elseif iscell(function_params{nd})
            this_function_params{np} = function_params{nd}{job_ndgrid{np}(nd)};
        elseif isstruct(function_params{nd})
            this_function_params{np} = function_params{nd}(job_ndgrid{np}(nd));
        end

    end
    this_function_params_string = 

    if iscompiled,
        
        %TODO :: IMPLEMENT COMPILED CODE
        
    else
        
        qsub_command_base = './matlab -nosplash -nodisplay -nojvm';
        qsub_command_file = sprintf('-o %s -e %s -N %s -r "%s(%s)"',...
            clusteringOptions.qsub_stdout_filename,      ...    % -o stdout file
            clusteringOptions.all_qsub_stderr_filename,  ...    % -e error file
            clusteringOptions.all_qsub_data_in_filename, ...    % -N job name (using data filename as proxy)
            function_to_call, ...    % data_in filename passed to subspace_clusterInstability
            this_function_params_string ...    % data_out filename passed to subspace_clusterInstability
            );
        qsub_command_resources = '-l num_cpus>15';
        qsub_command_full{i_oG} = [qsub_command_base ' ' qsub_command_file ' ' qsub_command_resources];
    end
    
end


%%%%%%%%%%%%%%%%%%%%%  END VAGUELY VALID CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


qsub_data_folder =  fullfile(analysis_save_directory,'qsub_data');
try mkdir(qsub_data_folder); catch,end

% Set up filenames for qsub data
all_qsub_data_in_filenames{i_oG}   = fullfile(qsub_data_folder , sprintf('qsub_data_in :: subspace_clusterInstablity :: %s :: %04.0f.mat',analysis_date,i_oG));
all_qsub_data_out_filenames{i_oG} = fullfile(qsub_data_folder , sprintf('qsub_data_out :: subspace_clusterInstablity :: %s :: %04.0f.mat',analysis_date,i_oG));
all_qsub_stdout_filenames{i_oG} = fullfile(qsub_data_folder , sprintf('qsub_stdout :: subspace_clusterInstablity :: %s :: %04.0f.out',analysis_date,i_oG));
all_qsub_stderr_filenames{i_oG}  = fullfile(qsub_data_folder , sprintf('qsub_error :: subspace_clusterInstablity :: %s :: %04.0f.err',analysis_date,i_oG));

% Save important filenames to the clusteringOptions struct
clusteringOptions.qsub_data_folder = qsub_data_folder;
clusteringOptions.qsub_data_in_filename = all_qsub_data_in_filenames{i_oG};
clusteringOptions.qsub_data_out_filename = all_qsub_data_out_filenames{i_oG};
clusteringOptions.qsub_stdout_filename = all_qsub_stdout_filenames{i_oG};
clusteringOptions.qsub_stderr_filename = all_qsub_stderr_filenames{i_oG};

% Generate and save all the inputs to
% subspace_clusterInstability for this run
qsub_data = { ...
    scoresForInstability,            ...
    clusteringOptions,               ...
    subsampleFraction,               ...
    randomRotationAmount, ...
    doShuffle,            ...
    noiseAmount,          ...
    angleNoise,           ...
    numIterations,                   ...
    instabilityCalculation            };
save(all_qsub_data_in_filenames{i_oG},qsub_data{:})

% Build the qsub command
% TODO :: fix inputs here - not sure what logfile is.
qsub_command_base = './matlab -nosplash -nodisplay -nojvm';
qsub_command_file = sprintf('-o %s -e %s -N %s -r "subspace_clusterInstability(''%s'',''%s'')"',...
    clusteringOptions.qsub_stdout_filename,      ...    % -o stdout file
    clusteringOptions.all_qsub_stderr_filename,  ...    % -e error file
    clusteringOptions.all_qsub_data_in_filename, ...    % -N job name (using data filename as proxy)
    clusteringOptions.all_qsub_data_in_filename, ...    % data_in filename passed to subspace_clusterInstability
    clusteringOptions.all_qsub_data_out_filename ...    % data_out filename passed to subspace_clusterInstability
    );
qsub_command_resources = '-l num_cpus>15';
qsub_command_full{i_oG} = [qsub_command_base ' ' qsub_command_file ' ' qsub_command_resources];

% Submit the job

fprintf('qsub_command_full{i_oG} :: \n %s \n',qsub_command_full{i_oG})

if do_qsub_submissions
    
    fprintf('SUBMITTING!\n')
    [status,cmdout{i_oG}] = system(qsub_command_full{i_oG});
    fprintf('Status :: %s\n\n',status)
    
else
    fprintf('Skipping qsub submission - run manually?\n')
    keyboard
end

fprintf('\n%s\n\tSubmitted%s\n\tOutput:\n%s\n',datestr(now,31),qsub_command_full{i_oG},cmdout{i_oG})