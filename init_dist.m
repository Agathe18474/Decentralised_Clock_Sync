function node_val = init_dist(type, n, min, max)
    % Initialse node values using given distribution
    % Inputs:
    %   type            type of distribution ("rdm" or "uniform" or "normal")
    %   n               nb of nodes
    %   min             minimum value
    %   max             maximum value
    % Output:
    %   node_val        array of values following given distribution

    range = max-min;
    
    switch type
        case "rdm"
            node_val = round(rand(n,1)*10)+min;
        case "uniform"
            node_val = min:(range+min);
        case "normal"
            mu = mean([min max]);
            r = max-mu; 
            node_val = random('Normal', mu, r/2, [n,1]);
    end
end