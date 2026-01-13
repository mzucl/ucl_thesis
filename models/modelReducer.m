% REDUCER FUNCTION
% Aggregates all emitted statistics from the mappers.
% Mixed reducer behavior:
%   - Sum the values for most keys
%   - Concatenate the values for specific keys (e.g., 'Z_E')
function modelReducer(~, intermValIter, outKV)
    % Initialize sumStats using the first value
    if hasnext(intermValIter)
        sumStats = getnext(intermValIter);
    else
        sumStats = struct(); % Return empty if nothing received
        return;
    end
    % Iterate through the rest of the intermediate values
    while hasnext(intermValIter)
        val = getnext(intermValIter);
        statFields = fieldnames(sumStats);

        for i = 1:numel(statFields)
            field = statFields{i};
            % if strcmp(field, 'E_Z')
            %     sumStats.(field) = [sumStats.(field), val.(field)];
            % else
            %     sumStats.(field) = sumStats.(field) + val.(field);
            % end
            sumStats.(field) = sumStats.(field) + val.(field);
        end
    end

    % Store the aggregated result
    add(outKV, 'aggregatedStats', sumStats);
end
