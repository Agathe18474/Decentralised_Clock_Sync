function [A, min_indeg_input, network_description] = connection_based_graphs(network_node_nb, k, conn_type)
    % Create random graph with specified number of random in-degrees 
    % INPUTS:
    % network_node_nb       nb of nodes in network
    % k                     minimum in-deg
    % conn_type             connectivity across network 
    %                           "same_conn" same connectivity at each node
    %                           "normal_dist_conn" conn at each node is a
    %                               normal distribution 
    % OUTPUTS:
    % A                     adjacency matrix
    % min_indeg_input       minimum in-degree
    % network description   [in_deg_ranks_min", in_deg_ranks_avg]
    %                           if "same_conn" used, the two vals will be
    %                           the same

max_indeg = network_node_nb-1;
min_indeg_input = k;    
n = network_node_nb;
A  = zeros(n, n);

% SAME CONNECTIVITY AT EACH NODE
if conn_type == "same_conn"
    for o=1:n
        sample_array = 1:n;
        sample_array(o) = [];
        rdm_sample =randsample(sample_array, min_indeg_input);
        A(rdm_sample,o) = 1;
    end

% CONNECTIVITY AT EACH NODE = NORMAL DIST BETWEEN MAX - MIN INDEG
elseif conn_type == "normal_dist_conn"
    % initialisation network parameters
    for o=1:n
        in_deg = min_indeg_input + round(rand*(max_indeg-min_indeg_input));
        %in_deg = min_indeg_input + floor(rand*(n-min_indeg_input));
        %index_rand = floor((n)*rand(in_deg,1)+1);
        sample_array = 1:n;
        sample_array(o) =[];
        rdm_sample =randsample(sample_array, in_deg);
        A(rdm_sample,o) = 1;
    end

    % if larger than in_deg input
    if nnz(sum(A,1)>min_indeg_input)==n
        rdm_array = randsample(1:n, 1);
        sample_array = 1:n;
        sample_array(rdm_array) =[];
        rdm_sample =randsample(sample_array, min_indeg_input);
        A(:, rdm_array) = 0;
        A(rdm_sample,rdm_array) = 1;
    end
% CONNECTIVITY AT EACH NODE = NORMAL DIST BETWEEN MAX - MIN INDEG
elseif conn_type == "random_conn"
    for o=1:n
        k = 3 + floor((n-3)*rand);
        sample_array = 1:n;
        sample_array(o) =[];
        sample_array = randsample(sample_array, k);
        A(sample_array,o) = 1;
    end
end


G = digraph(A);
H = simplify(G);

in_deg_ranks = centrality(H,'indegree');
in_deg_ranks_avg = mean(in_deg_ranks);
in_deg_ranks_min = min(in_deg_ranks);

A = adjacency(H);

network_description = [in_deg_ranks_min, in_deg_ranks_avg];

end