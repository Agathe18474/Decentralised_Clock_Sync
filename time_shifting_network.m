clearvars 
close all 

% Simulation comparing filter types on single run 
%
% DESCRIPTION:
%   Comparing filter types for specific disruptive behaviours on random or
%   known graph scenarios. 
%
% STRENGTH
%   Fault resilience (up to UNKNOWN faults) 
%   No Network topology knowledge is required
%   Knowledge of number of faults in system not required 
%
% WEAKNESSES
%   Synchronous/no delay tolerance
%   No dynamic byzantine behaviour (no colluding)
%   Non intelligent byzantine agents
%   No time shifting tolopogy 
%   No drift 
%   POSSIBLE: clustering might not work for low values
%
% CALLED FUNCTIONS:
%   disruptor_parameters
%   disruptor_node_val
%   init_dist
%   connection_based_graphs
%   graphs_test_cases
%   update_script
%   opinion_plots_single
%
% AUTHOR:
%  Agathe BOUIS 12/02/2024
%
% INPUTS
% experiment_1              run experiments from Bouis, A., Clark R. A., 
%                               Macdonald M, (2024) "Engineering consensus 
%                               in static networks with unknown disruptors" 
%                                   > "1a" - experiement 1a
%                                   > "1a" - experiement 1b
%                                   > "null" - select inputs for simulation 
%                                       of choice                          
% filter_type_array         type of filter tested
%                               "ODDI-C"
%                               "MSR"
%                               "Mean-based"
% rdm_graph_flag            is the graph tested randomly created? [boolean]
%                               > true: rdm graph. if true, specify:
%                                   max_conn, conn_type, min_indeg_input
%                               > false: specify graph test case
%                                   name
% create a random graph
% n                         number of nodes in network tested
%       max_conn                    max connectivity (in-deg) tested
%       conn_type                    type of graph tested 
%                                       "same_conn" - same connectivity for 
%                                       each node
%                                       "normal_dist_conn" connectivity based 
%                                       on normal distribution
% use already created graph
%       min_indeg_input             min in-deg allowed (if not same conn)
%       name                        test scenario investigated
%                                       "Test_1": 7 nodes (2,2)-robust graph
%                                       "Test_2": 8 nodes (3,1)-robust graph
%                                       "Test_3": 14 nodes (2,2)-robust graph
%                                       "Test_4": 7 nodes (3,3)-robust graph 
%                                       "Test_5": 5 nodes (2,2)-robust DIgraph
%                                       "Test_6": 6 nodes (3,1)-robust DIgraph
%                                       "Test_7": 12 nodes (3,2)-robust DIgraph
%                                       "Test_8": 15 nodes (5,1) or (4/3,2)-robust DIgraph
%                                       "Test_9": 7 nodes (4,1) or (3,2)-robust DIgraph                

% initial node value distribution 
% init_dist_type            initial distribution type
%                               "rdm": random dist
%                               "uniform": uniform 
%                               "normal": normal 
% max_allowable             maximum allowable in distribution 
% min_allowable             min allowable in distribution 

% simulation parmeters: 
% max_time_step             max possible nb of time steps
% break_time                what run scenario [boolean]
%                               true: let run until stable or meet max time
%                               false: run until max time
% mu                        learning rate/update step
% abs_zero_val              min abstolute difference
%                               calculate min value of convergence metric 
%                               from it
%                               replace "0" by - usually 10^(-7)

% plot_comp_iteration       choice of which plots to output [boolean]
%                               true/false - "Node Trajectory" 
%                               true/false - "Convergence Metric"
%
%
% OUTPUTS
% node_val                  evolution of node values over time 
% t                         associated time array
% convergence metric        relative convergence between compliant nodes
% POSSIBLE PLOTS
% Per each connectivity test [plot_comp_iteration]
%   node opinion trajectory
%   convergence metric

%% INPUTS 
% load and use satellite data?
% if satellite data not used, random network with rolling adjacency matrix
% created
sat_data_tf=false;
sat_data_name= "ISL_Contacts_10shells500sats.mat";

n = 20;
n_temp = 20;%round(n/3);
%n = 500;

% disruptive behaviour
% either manual input, or provide the number of disruptors wanted
% manual 
 i_dis = [];
% random generation
% i_dis_nb = 50;
% sample_array = 1:n;
% i_dis =randsample(sample_array, i_dis_nb);
% i_dis = sort(i_dis);

% node values initial 
init_dist_type = "normal";
min_allowable = 1;
max_allowable = 20;
% simulation parameters
max_time_step = 500; 
break_time = false;         
mu = 0.5; 
test_step_size = 1;
abs_zero_val = 10^(-7);
% plotting directory
plot_dict_comp_iteration = [true false];
% filter type - COMPARISON
%filter_type_array = ["MSR", "ODDI-C"];
filter_type_array = ["ODDI-C"];

%% INITIALISATION 
% Actual length of test
comp_test = filter_type_array;
possible_dis_types = ["sine_wave","lin_const", "noise"];

% if use satellite data
if sat_data_tf 
    sat_data=load("ISL_Contacts_10shells500sats.mat");
    A = sat_data.ISL_Contacts;
    n = size(cell2mat(A(1)),1);
else
    A = rolling_adjacency_matrix(n,n_temp);
end


dis_type = randsample(possible_dis_types, length(i_dis), true);
dis_param = disruptor_parameters(dis_type, max_allowable, min_allowable, false, max_time_step);

dis_val=[];
for dis_id=1:length(i_dis)
    dis_val(dis_id,:) =  disruptor_node_val(dis_type(dis_id), dis_param(dis_id, :), max_allowable, min_allowable, (1:max_time_step+2));
end

convergence_metric_comparison = zeros(length(comp_test), max_time_step+2);

% start values
node_val = init_dist(init_dist_type, n, min_allowable, max_allowable);
zero_val = abs_zero_val/(max(node_val)-min(node_val));
     
for comparison_test=1:length(filter_type_array)
    filter_type = filter_type_array(comparison_test);
    percentage_filtered_comp=zeros(n-length(i_dis), 1);
 
    % update 
     [t, node_val, diff_nb, percentage_filtered_comp] = update_script_time_shifting(filter_type, max_time_step, node_val(:,1), ...
     A, mu, i_dis, dis_val, break_time, percentage_filtered_comp);
    
    % calculate & store convergence metric
    convergence_metric_comp_iteration = diff_nb(1,:)./max(diff_nb(1,:));
    convergence_metric_comp_iteration(convergence_metric_comp_iteration<zero_val)=zero_val;
    convergence_metric_comparison(comparison_test, :) = convergence_metric_comp_iteration;
    
    % plot single run
    opinion_plots_single(node_val, t, i_dis, comparison_test, ...
        plot_dict_comp_iteration, convergence_metric_comp_iteration, zero_val)
end

    
f=figure();
t_cm = repmat(t,length(comp_test),1);
semilogy(t_cm.',convergence_metric_comparison.') 


hold on
ax = gca; 
styles = ["-","--",":", "-."];
ax.LineStyleOrder = styles; 
xlim([0 t(end)])
ylim([zero_val 1]) 
grid on
legend(filter_type_array)
xlabel("Time steps")
ylabel("Convergence Metric")
grid on
box on 
fontsize(f,13, "points")