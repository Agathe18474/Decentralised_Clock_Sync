function break_check = loop_check(node_values_new, node_values_current, node_values_past, error_thresh)
    % Check if continue looping of if values have converged normalised node difference
    % Inputs:
    %   node_values_new         new value of nodes
    %   node_values_current     current value of nodes
    %   node_values_pas         past value of nodes
    % Output:
    %   break_check             boolean, if values stable=1, else=0

    % check if the node values have converge to 1 value
    % if the last two iterations have the same val as the current value, 
    % break out of loop
    current_round = abs(round(node_values_new,2) - round(node_values_current,2)) <= error_thresh;
    next_round = abs(round(node_values_current,2)-round(node_values_past,2)) <= error_thresh;

    if all(current_round) && all(next_round)
        break_check = true; 
    else 
        break_check = false;
    end 
end
