%initialise_ASYNC
close all 
clear all

% simulation parameters 
max_time_step = 500; 
break_time = false;  
alpha = 0.8;
abs_zero_val = 10^(-9);
conn_type = "Sat Network";
%conn_type = "Sat Network Josh";
% conn_type = "Rings";
nb_rings = 2;
tof_zeroing = true;
delay_correction = false;
synchronous_check = false;
event_trigger_check = false;
noise_check = false; 

manual_input_values = false;

% monte-carlo length 
monte_carlo = 2;
i_dis_nb = 0;

% filter type 
% filter_type = "ODDI-C";
% filter_type = "Hampel";
% filter_type = "Dixon";
filter_type = "Mean-based"; 

% cycles (Consensus, Broadcast, Topology Change, Prune)
CC = 10;
BC = 3;
TC = 3;
PC = CC;

% works for network sat
% gamma = 0.8;
% eta = gamma^2;
gamma = 1;
eta = 0.5;

if tof_zeroing
    delay_correction = false;
end

% based on input connectivity type - load data 
if conn_type == "Sat Network Josh"
    data = "ISL_Contacts_10shells200sats.mat";
    sat_data=load(data);
    A = sat_data.ISL_Contacts;
    A1 = A(1, 100:end);
    n = size(cell2mat(A(1)),1);
    tof = 0.0016*ones(n);
    conn_start= 1; 
    conn_max=1;
    test_step_size = 1;
    conn_test = conn_start:test_step_size:conn_max;
    names = "";
    % set rings to 0 if forgotten
    nb_rings = 0;
elseif conn_type == "Sat Network"
    data = "200_sats_random_network.mat";
    %data = "500_sats_random_network.mat";
    sat_data=load(data);
    A = sat_data.A;
    tof = sat_data.tof;
    n = sqrt(size(A,1));
    conn_start= 1; 
    conn_max=1;
    test_step_size = 1;
    conn_test = conn_start:test_step_size:conn_max;
    names = "";
    % set rings to 0 if forgotten
    nb_rings = 0;
elseif conn_type == "Rings"
    [A, tof, n, sat_list_names, names, conn_start,test_step_size, conn_max] = ring_setup(nb_rings);
    conn_test = conn_start:test_step_size:conn_max;    
end

% If want to investigate a SINGLE ring 

% data = "2_rings_6_conn_ANALYSIS.mat";
% sat_data=load(data);
% A = sat_data.A;
% tof = sat_data.tof;
% pos_plot = sat_data.pos_plot;
% n = sqrt(size(A,1));
% %tof = 0016*ones(n);
% conn_start= 2; 
% test_step_size=2;
% conn_max = 10;
% conn_test = conn_start:test_step_size:conn_max;
% A_type = class(A);
% A = time_step_adajacency(A, A_type,1, n, tof);
% A = double(A);

% disruptive behaviour
% either manual input, or provide the number of disruptors wanted
% manual 
% i_dis = [];
% random generation
sample_array = 1:n;
i_dis = randsample(sample_array, i_dis_nb);
i_dis = sort(i_dis);

% node values initial 
init_dist_type = "normal";
min_allowable = 1;
max_allowable = 3;

% plotting directory
% plot_dict_single - what to plot for 1 monte-carlo
% plot_dict_multi - what to plot when looking at avg of mc 
plot_dict_single = [true false false];
plot_dict_multi = [false true false];
plot_colour = "multi";

colour_id = ["#0072BD", "#D95319", "#EDB120", "#7E2F8E", "#77AC30"];
style_id = {'-','--', ':', '-.'};

possible_dis_types = ["sine_wave","lin_const", "noise_const"];

% if too many nodes, force monochromatic colouring
if n>40
    plot_colour = "mono";
end


% if using manual input
% simulation parameter initialisation
% [clock_parameters, node_val, zero_val] = simulation_init(n, init_dist_type, min_allowable, max_allowable, abs_zero_val);
% simulation_data = readmatrix("200_node_network.xlsx");
% clock_parameters = simulation_data(:,1:4);
% node_val = simulation_data(:,end);

% initialise arrays 
convergence_metric_rate_mc = zeros(monte_carlo, max_time_step+1);
convergence_metric_rate = zeros(length(conn_test), max_time_step+1);
convergence_metric_rate_max = zeros(length(conn_test), max_time_step+1);
convergence_metric_rate_min = zeros(length(conn_test), max_time_step+1);

convergence_metric_time_mc = zeros(monte_carlo, max_time_step+1);
convergence_metric_time = zeros(length(conn_test), max_time_step+1);
convergence_metric_time_max = zeros(length(conn_test), max_time_step+1);
convergence_metric_time_min = zeros(length(conn_test), max_time_step+1);

counter = 1;

tic
for test=1:1
    for conn_iteration=conn_start:test_step_size:conn_max
        for mc_iteration=1:monte_carlo 
            % change initial node conditions for each run
            [clock_parameters, node_val, zero_val] = simulation_init(n, init_dist_type, min_allowable, max_allowable, abs_zero_val);  
 

            if tof_zeroing==true
                tof = zeros(n);
            end
            % if use manual input - override 
            if manual_input_values
                init_data=load("200_sats_init_2.mat");
                node_val = init_data.node_val_init;
                clock_parameters = init_data.clock_param_init;
            end
            
            % disruptive nodes creation
            dis_type = randsample(possible_dis_types, length(i_dis), true);
            dis_param = clock_disruptor_parameters(dis_type, max_allowable, min_allowable, false, max_time_step);
            dis_val=zeros(length(i_dis), max_time_step+1);
            for dis_id=1:length(i_dis)
                dis_val(dis_id,:) =  clock_disruptor_val(dis_type(dis_id), dis_param(dis_id, :), max_allowable, min_allowable, (1:max_time_step+1));
            end

          
            % simulation
             [t, node_val, absolute_rate, relative_time_diff, relative_rate_diff, max_error_time, max_error_rate, delay] = script_ASYNC(filter_type, max_time_step, ...
            node_val, A, i_dis, dis_val, clock_parameters, alpha, tof, CC, BC, TC, PC, eta, gamma, synchronous_check, delay_correction, event_trigger_check, noise_check);
           

            % calculate & store convergence metric - per monte carlo iteration
            % ie, simulating the same graph (same conn) multiple times          
            % convergence in time & rate
            convergence_metric_time_mc(mc_iteration, :) = convergence_metric_single(relative_time_diff, zero_val);
            convergence_metric_rate_mc(mc_iteration, :) = convergence_metric_single(relative_rate_diff, zero_val);

            
            % if first mc iteration, plot simulation output
            if mc_iteration == 1
                % Time plots
                opinion_plots_single(node_val, t,convergence_metric_time_mc, ...
                    max_error_time, i_dis, plot_dict_single, ...
                    zero_val,plot_colour, "Time") 

                % Rate plots
                opinion_plots_single(absolute_rate, t,convergence_metric_rate_mc, ...
                    max_error_rate, i_dis,plot_dict_single, ...
                    zero_val,plot_colour, "Rate") 
            end

        end
    
        % monte carlo analysis 
        % calculate & store convergence metric - average metric across the
        % mc of the same graph (same conn)
        % in time & rate
        [convergence_metric_time(counter,:),convergence_metric_time_min(counter,:),convergence_metric_time_max(counter,:)]  = convergence_metric_multi(convergence_metric_time_mc);
        [convergence_metric_rate(counter,:),convergence_metric_rate_min(counter,:),convergence_metric_rate_max(counter,:)]  = convergence_metric_multi(convergence_metric_rate_mc);            
           
            
        if conn_iteration<conn_max        
            % re initialise initial node values, zero-val, & clock param with 
            % new n (for ring nb>1, n increases w/ connectivity) 
            if nb_rings==1
                altitude = 2000;
                rE = 6371;
                r = (altitude+rE)*10^(3);
                k = conn_test(counter);
                [A, tof] = ring_graph(r, n, k);
            elseif nb_rings>1
                % load new sat data 
                data = sat_list_names(counter);
                sat_data=load(data);
                A = sat_data.A;
                tof = sat_data.tof;
                % A = sat_data.A_temp;
                % n = sqrt(size(A,1));
                % tof = zeros(n,length(A));
                
            end
            % force tof to zero - instant comms
            if tof_zeroing==true
                tof = zeros(n);
            end
            counter = counter+1; 
            %[clock_parameters, node_val, zero_val] = simulation_init(n, init_dist_type, min_allowable, max_allowable, abs_zero_val);
        end
    end
    
    % create plot names
    if isempty(names)
        name_spec = "Sat n= %G";
        for i=1:length(conn_test)
            names(1,i) = sprintf(name_spec,n);
        end
    end    

    boxplot_data = [];
    percentage_filtered= [];
    percentage_filtered_max=[];
    percentage_filtered_min = [];

    % PLOT MULTIPLE RUNS
    if sum(plot_dict_multi)~=0

        if length(conn_test)==1
            names = {"Time", "Rate"};   
            convergence_metrics = [convergence_metric_time;convergence_metric_rate];
            convergence_metrics_max = [convergence_metric_time_max;convergence_metric_rate_max];
            convergence_metrics_min = [convergence_metric_time_min;convergence_metric_rate_min];


            % Time and rate plotted together
            opinion_plots_multi_combined(t, convergence_metrics, ...
            convergence_metrics_max, convergence_metrics_min,  names)
        else
            % Time
            opinion_plots_multi(t, convergence_metric_time, ...
                 convergence_metric_time_max, convergence_metric_time_min, boxplot_data, percentage_filtered, ....
                 percentage_filtered_max, percentage_filtered_min, plot_dict_multi, ...
                 filter_type, names, zero_val, "Time")
    
             % Rate
             opinion_plots_multi(t, convergence_metric_rate, ...
                 convergence_metric_rate_max, convergence_metric_rate_min, boxplot_data, percentage_filtered, ....
                 percentage_filtered_max, percentage_filtered_min, plot_dict_multi, ...
                 filter_type, names, zero_val, "Rate")
        end
    end
    toc
end