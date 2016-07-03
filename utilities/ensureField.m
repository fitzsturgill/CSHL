function s = ensureField(s, f, type)
% add field to a nonscalar structure array of a specified type
% s- structure f- fieldName type: 'mat', 'cell', 'struct', 'NaN'
% if field preexists, nothing is done
    if nargin < 3
        type = 'mat';
    end
    
    if ~isfield(s, f)
        [s(:).(f)] = deal(cell(size(s)));
        for counter = 1:numel(s)
            switch type
                case 'mat'
                    s(counter).(f) = [];
                case 'cell'
                    s(counter).(f) = cell();
                case 'struct'
                    s(counter).(f) = struct();
                otherwise
            end
        end
    end