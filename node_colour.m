function node_colour_idx = node_colour(node_value, max_init, min_init, nb_colours)
    % Node colours based on its current value wrt to initial t(0) values   
    % Inputs:
    %   node_value          value of nodes at given time step
    %   min_init            minimum initial value
    %   max_init            maximum initial value
    % Output:
    %   node_colour_idx     colours of each node

    node_colour_idx = round(((node_value-min_init)/(max_init-min_init))*(nb_colours-1))+1;
end