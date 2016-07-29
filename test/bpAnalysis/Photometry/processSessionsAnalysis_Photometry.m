function sessions = processSessionsAnalysis_Photometry(sessions, analysisName, varargin)
    if nargin < 2 || isempty(analysisName)
        analysisName = 'analysis'; % default behavior, see processAnalysis_Photometry
    end
    
    if ischar(analysisName)
        if ~isfield(sessions, analysisName)
            [sessions(:).(analysisName)] = deal(cell(size(sessions))); % use ensureField in future for this...
            for counter = 1:length(sessions)
                sessions(counter).(analysisName) = struct();
            end
        end
    end
    
    if nargin < 3
        varargin = []; % again default behavior
    end

    for counter = 1:length(sessions)
        disp(['*** Processing: ' sessions(counter).filename '***']);       
        sessions(counter).(analysisName).Photometry = processAnalysis_Photometry(sessions(counter), varargin{:});
    end