function max_difference = max_difference(data, i_non_dis)
    % maximum error between nodes for each time step 
    % Inputs: 
    %   data                    matrix of data (time or rate) for each node
    %                               at each time step
    % Outputs: 
    %   max_difference          calculated maximum difference  

    max_val = max(data(i_non_dis,:));
    min_val = min(data(i_non_dis, :));

    max_difference = max_val-min_val;

end