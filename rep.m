function x = rep(x, oldvalues, newvalues)
% Function to replace numbers in x (oldvalues) by new numbers (newvalues).
% If x is a cell, apply the replace routine to each element of the cell.
if iscell(x)
    N_cell_elements = numel(x);
    for n = 1:N_cell_elements
        x{n} = rep_routine(x{n}, oldvalues, newvalues);
    end
else
    x = rep_routine(x, oldvalues, newvalues);
end
end

function x = rep_routine(x, oldvalues, newvalues)
if ~isempty(x)
    [~, locs] = ismember(x, oldvalues);
    nonzerolocs = locs ~= 0;
    x(nonzerolocs) = newvalues(locs(nonzerolocs));
end
end