function diff_nb = node_diff(node_values, max_init, min_init)
    % Calculate normalised node difference
    % Inputs:
    %   node_values         array of values from non byzantine nodes
    %   max_init            maximum initial value
    %   min_init            minimum initial value
    % Output:
    %   diff_nb             normalised difference between the nodes

    % convergence values in round
    convergence = round(node_values,2);
    % difference between the nodes
    diff = round((convergence - convergence.'),3);
    diff = triu(diff);
    diff_norm = diff/(max_init-min_init);
    diff_norm = diff_norm.^2;
    diff_rows = sum(abs(diff_norm));
    nb_diff = 0;
    for i = 1:length(diff)-1
        nb_diff = nb_diff + i;
    end
    diff_nb = sum(diff_rows)/nb_diff; 

end