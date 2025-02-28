function [clock_parameters, node_val, zero_val] = simulation_init(n, init_dist_type, min_allowable, max_allowable, abs_zero_val)
    % Initialisation of simulation parameters
    % Clock parameters, Node Init V-times, Relative zero val
    % Inputs: 
    %   n                   # nodes
    %   init_dis_type       distribution of nodes' initial values
    %   min_allowable       min value allowed for initial value
    %   max_allowable       max value allowed for initial value
    %   abs_zero_val        absolute zero - value past which convergence
    %                            metric is considered zero 
    % Outputs: 
    %   clock_parameters    clock parameters
    %   node_val            node values
    %   zero_val            zero value (floor) relative to initial node
    %                           values
    
    clock_parameters = clock_param(n); 
    node_val = init_dist(init_dist_type, n, min_allowable, max_allowable);
    zero_val = abs_zero_val/(max(node_val)-min(node_val));

end