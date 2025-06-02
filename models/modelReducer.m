% REDUCER (sum all emitted stats)
% Mixed reducer: some keys you want to sum, and one key (e.g. 'Z_E') you want to concatenate
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
            % Concatenate matrices horizontally
            if strcmp(field, 'E_Z')
                sumStats.(field) = [sumStats.(field), val.(field)];
            else
                sumStats.(field) = sumStats.(field) + val.(field);
            end
        end
    end

    % Store the aggregated result
    add(outKV, 'aggregatedStats', sumStats);
end
