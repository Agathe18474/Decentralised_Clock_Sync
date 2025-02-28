function [t, V, absolute_rate ,relative_time_diff, relative_rate_diff, max_error_time, max_error_rate,stored_delay_confirmed] = script_ASYNC(filter_type, max_time, ...
    node_val_init, A, i_dis, dis_val, clock_parameters, alpha, tof, CC, BC, TC, PC, eta, gamma, synchronous_check, delay_correction, event_trigger_check, noise_check)
    % Asynchronous Update Script 
    %
    % DESCRIPTION: 
    %       
    % AUTHOR:
    %   Agathe BOUIS, 05/09/2024
    %
    % INPUTS:
    %   filter_type         type of filter used 
    %                           > ODDI-C
    %                           > MSR
    %                           > Hampel
    %                           > Dixon
    %                           > Mean-based
    %   max_time            maximum nb of time steps before break 
    %   node_val_init       value of nodes at time 0 (initialised value)
    %   A                   3D adjacency matrix of graph (full matrix)
    %                       3rd dimension is time variation
    %   i_dis               id of disruptive nodes
    %   dis_val             value of disruptive nodes for each time step
    %   clock_parameters    H-clock parameters
    %   alpha               smoothing factor
    %   tof                 time of flight between nodes for each topology
    %   CC                  Consensus Cycle: how often node updates its rate
    %                           given as constant
    %                           calculated as multiple of dH (H-tick)
    %   BC                  Broadcast Cycle: how often node broadcasts
    %                           given as constant
    %                           calculated as multiple of dH (H-tick)
    %   PC                  Prune Cycle: how long data is kept for update
    %                           given as constant
    %                           calculated as multiple of dH (H-tick)
    %  TC                   Topology Change: nb time steps before change
    %                           topology (next Adjacency matrix)
    %
    % OUTPUTS:           
    %   t                   time elapsed
    %   V                   array of node V-times for each time step
    %   relative_time_diff  normalised relative time difference
    %   relative_rate_diff  normalised relative rate difference
    %   max_error_time      max absolute error in time
    %   max_error_rate      max absolute error in rate
    %
    % CALLED SCRIPTS:
    %   node_diff  
    %   filter_ODDI_C
    %   filter_MSR
    %   loop_check

    % ADD CDI to code
    % path creation to use Communities of Dynamical Influence scripts
    % Determine where your m-file's folder is.
    % folder = fileparts(which(mfilename)); 
    % % Add that folder plus all subfolders to the path.
    % addpath(genpath(folder));
    
    % format initial V-time
    V = node_val_init; 

    % safe guards
    smallest_rate=0.01;
    max_size_database = 300;
    %max_size_database = 100;
    min_size_database = 10;
    min_nb_sources = 4;


    % clock parameters
    R=clock_parameters(:,1);
    O=clock_parameters(:,2);
    freq=clock_parameters(:,3);     
    
    % how are the two hop data points sampled?
    database_2_hop_sampling = "random";
    % database_2_hop_sampling = "median";
    % can also sample "median" for each ID
    
    % initialise values for start of loop
    break_check = false;
    k = 1;    
    t = 1;
    n = length(V);
    i_normal = 1:n;
    if ~isempty(i_dis)
        i_normal(i_dis) = [];
    end
    
    % initialise arrays
    rate = zeros(n,1);
    rate(i_normal) = ones(length(i_normal),1);    
    last_update = zeros(n,1);
    size_data_base = zeros(n,1);
    delay = zeros(n,1); 
    delay_rate = ones(n,1);
    delay_confirmed = zeros(n,1);
    tolerance = 0.1; 
    overlap = ones(n,1);
    outlier_check = false(n,1);
    V_last_broadcast = inf(n,1);
    percentage_filtered = nan(n,1);
    
    C = zeros(n,1);
    D = zeros(n,1);
    C_rate = rate(:,1);
    t_TC = 2;
    TC_counter = 1;
    Database_Stored = cell(n,1);
    Database_Received = cell(n,1);
    Database_New = cell(n,1);    

    % initialise databases
    for i=1:n  
         Database_Stored{i}=double.empty([0 6]);
         Database_Received{i}=double.empty([0 6]);
         Database_New{i}=double.empty([0 6]);
    end


    % initialise H-times & triggers
    dH = H_clock(R, O, freq);  
    H = zeros(n,1);
    dH = H_clock(R,O,freq);    
    
    % if synchronous, force BC=CC
    if synchronous_check
        BC = CC;
        H_update = dH.*CC.*ones(n,1);
        %H_update = CC.*ones(n,1);
        H_update(i_dis) = inf;
        H_broadcast = zeros(n,1);
    else
        H_update = 3.*BC.*dH.*rand(n,1);
        %H_update = 2.*BC.*rand(n,1);
        %H_update=CC+H;
        H_update(i_dis) = inf;
        %H_broadcast = BC.*dH.*rand(n,1);
        H_broadcast = BC.*rand(n,1);
        %H_broadcast = BC+H;
    end
    
    % initial nodes' H- & V-times
    absolute_rate(:,k) = rate(:,k).*dH;

    % find min and max values among nodes
    max_init_time = max(V(:,1));
    min_init_time = min(V(:,1));
    max_init_rate = max(absolute_rate(:,1));
    min_init_rate = min(absolute_rate(:,1));

    % adapted update_parameter
    CC_adapted = CC*ones(n,1);

    % intial relative differences in time & rate between nodes
    relative_time_diff = node_diff(V(i_normal,k), max_init_time, min_init_time);
    relative_rate_diff = node_diff(V(i_normal,k), max_init_rate, min_init_rate);
   
    A_type = class(A);
    
    % in-test: if ring, min nb sources should be ring conn + 1
    if contains(A_type, 'double')
        min_nb_sources = mode(sum(A,1))+1;
    else
        A_current_time = time_step_adajacency(A, A_type, 1, n, tof);
        %min_nb_sources = mode(sum(A_current_time, 1))+1;        
    end

    
    % SIMULATION
    while ~break_check
        % increase time counters
        k = k+1;
        % universal/real time elapsed
        % increase time counters        
        t = cat(2,t,k);
        dH = H_clock(R,O,freq);
        H(:, k) = H(:,k-1)+dH;

        percentage_filtered = cat(2, percentage_filtered, nan(n,1));

        

        % if include noise in data
        if noise_check
            % power in db - 0.25 W = -6dB
            dw = wgn(n,1,-13);
            V(i_normal,k) = V(i_normal,k-1) + rate(i_normal,k-1).*dH(i_normal) + dw;
        else
            V(i_normal,k) = V(i_normal,k-1) + rate(i_normal,k-1).*dH(i_normal);
        end
        rate(:,k) = rate(:,k-1);
        absolute_rate(:,k) = rate(:,k-1).*dH;
        constant_delay_1_hop(:,k)=delay;
           
        
        % if Topology Change trigger time 
        if k-t_TC>TC || k==2
            [A_current_time, tof_current_time] = time_step_adajacency(A, A_type, TC_counter, n, tof);
            indegs(k,:) = sum(A_current_time, 1);
            t_TC = k;
            TC_counter = TC_counter +1;
        end
        
        % if disruptive nodes, disruptive nodes update based on
        % pre-calculated value
        if i_dis ~=0
            V(i_dis,k) = dis_val(1:length(i_dis),k); 
        end
                
        % check for BROADCAST triggers        
        if event_trigger_check
            broadcasting_nodes = event_triggering(V(:, k), V_last_broadcast);
        else
            broadcasting_nodes = H(:,k)>=H_broadcast;
        end
        id_broadcasting = find(broadcasting_nodes);

        for i=1:length(id_broadcasting)
            % broadcasting node ID
            node_id = id_broadcasting(i);
            
            % data to broadcast
            [data_broadcast, ID_broadcast] = DQF(Database_Stored{node_id}, V(node_id, k), "broadcast", delay_confirmed(node_id));
            
            % indexes of connected nodes
            conn_idx = find(A_current_time(node_id,:)).';
            
            % for each connected node, send data to their Receiving
            % database
            for j=1:length(conn_idx)
                % receptor node ID   
                node_id_receptor = conn_idx(j);
                
                % tof in terms of the H- & V-times of the receptor
                % [tof_H_receptor, tof_V_receptor]
                %% VERIFY
                tofs_receptor = tof_relative(tof_current_time(node_id, node_id_receptor), H(node_id_receptor, k)-H(node_id_receptor,k-1), V(node_id_receptor, k)-V(node_id_receptor,k-1));
                
                %% VERIFY
                % V_receptor = V(k, node_id_receptor);
                % tof_delay = tof_current_time(node_id, node_id_receptor);
                % formatting of received data by receptor node
                %Database_Appendage = database_appender(data_broadcast, ID_broadcast, V(node_id_receptor, k), H(node_id_receptor, k), tofs_receptor, V(node_id, k), node_id);               
                Database_Appendage = database_appender(data_broadcast, ID_broadcast, V(node_id_receptor, k), H(node_id_receptor, k), tofs_receptor, V(node_id, k), node_id, outlier_check(node_id));
                % append new data to receptor node's database
                Database_Received{node_id_receptor} =  [Database_Received{node_id_receptor}; Database_Appendage];                  
            end
            
            % update the broadcast trigger 
            if event_trigger_check
                V_last_broadcast(node_id) = V(node_id,k);
            else
                H_broadcast(node_id) = H_trigger(H(node_id,k), BC, dH(node_id), "future");
            end
        end


        % check for UPDATE triggers
        %database_size_check = and(Database_Received_size>max_size_database, delta_H<CC);
        %updating_nodes = or(H(:,k)>=H_update,database_size_check);
        Database_Received_nb_sources = cellfun(@(x) sum(x(:,5)), Database_Received);
        Database_Received_size = cellfun(@length, Database_Received);
        data_base_received(:,k) = Database_Received_size(:);

        [updating_nodes, change_update_pattern] = update_trigger(H(:,k), H_update, Database_Received_size, Database_Received_nb_sources, min_size_database, max_size_database, min_nb_sources);
        id_updating = find(updating_nodes);  
        
        % H-time since last update
        delta_H = H(:,k)-last_update;
        
        indices_max = 1:n;
        node_with_max_val = indices_max(V(:,k)==max(V(:,k)));

        for i=1:length(id_updating)
            % updating node ID
            node_id = id_updating(i);
            
            % Break in code in case of error 
            if k>200 && node_id==node_with_max_val % any(absolute_rate(:,k)>3*std(absolute_rate(:,k)))%
                oo=1;
            end  

            if ~isempty(Database_Stored{node_id})
                % before pruning, account for all IDs 
                ID_past = Database_Stored{node_id}(:,2);

                % prune older data
                % assume that PC is CC
                H_prune = H_trigger(H(node_id,k), delta_H(node_id), dH(node_id), "past");                
                %H_prune = H_trigger(H(node_id,k), CC_adapted(node_id), dH(node_id), "past");                
                Database_Stored{node_id} = database_prune(Database_Stored{node_id}, H_prune);
            end

            % update the stored database to contain: 
            %   > older, still relevant data (if PC>CC)
            %   > newer, filtered data
            Database_New{node_id} = [Database_Stored{node_id}; Database_Received{node_id}];

            % Database management 
            % filter database of received values
            % keep 1-hop transmission over 2-hop transmission from same node
            % keep 2-hop transmission received from node (either rdm or 
            % median for each ID)              
            % [Database_Filtered, nb_occurences_per_ID, nb_occurences_per_2_hop_ID] = database_filter(Database_Received{node_id}, node_id, "Step 1", database_2_hop_sampling);
            [Database_Stored{node_id}, nb_occurences_per_ID] = database_filter(Database_New{node_id}, node_id, "Step 1", database_2_hop_sampling);
            
            % empty reception database 
            Database_Received{node_id} = double.empty([0 6]);
            Database_New{node_id} = double.empty([0 6]);
            
            % CHANGE FROM PREVIOUS VERSION TO SPEED UP CODE
            % if ~isempty(Database_Stored{node_id})
            %     % before pruning, account for all IDs 
            %     ID_past = Database_Stored{node_id}(:,2);
            % 
            %     % prune older data
            %     % assume that PC is CC
            %     %H_prune = H_trigger(H(node_id,k), delta_H(node_id), dH(node_id), "past");                
            %     H_prune = H_trigger(H(node_id,k), CC_adapted(node_id), dH(node_id), "past");                
            %     Database_Stored{node_id} = database_prune(Database_Stored{node_id}, H_prune);
            % end
            % 
            % % update the stored database to contain: 
            % %   > older, still relevant data (if PC>CC)
            % %   > newer, filtered data
            % Database_Stored{node_id} = [Database_Stored{node_id}; Database_Filtered];
            % 
            % if size(Database_Stored{node_id},1)~=size(Database_Filtered,1)
            %     % filter database again to only keep latest value
            %     % if received any ID again, remove older data
            %     [Database_Stored{node_id}, nb_occurences_per_ID] = database_filter(Database_Stored{node_id}, node_id, "Step 2", database_2_hop_sampling , nb_occurences_per_2_hop_ID);
            %     %nb_occurences_per_ID = [nb_occurences_per_1_hop_ID;nb_occurences_per_2_hop_ID];
            % end

            size_data_base(node_id) = size(Database_Stored{node_id},1);

            %Database_Stored{node_id} = database_filter(Database_Stored{node_id});

            % consensus has anydata to pull from that database is empty
            if isempty(Database_Stored{node_id})
                database_check = true;
            else
                database_check = false;
            end
            
            
            if exist('ID_past','var')
                % account for all IDs 
                ID_current = Database_Stored{node_id}(:,2);
                % check overlap between current & past ID
                overlap_new =  ID_overlap(ID_past,ID_current);
                overlap(node_id) = single_exponential_smoothing(alpha, overlap_new, overlap(node_id));
            else
                overlap_new = 1;
                overlap(node_id) = single_exponential_smoothing(alpha, overlap_new, overlap(node_id));
            end
            
            
            % if correct for delays
            if delay_correction
                % calculate delays
                [delay(node_id), delay_rate(node_id), check_delay] = delay_calculation(Database_Stored{node_id}, delta_H(node_id), delay(node_id), alpha, delay_rate(node_id), min_nb_sources);
            
                % calculated delay correction from constant error 
                delay_confirmed(node_id) = delay_confirmation(delay_confirmed(node_id), delay(node_id), delay_rate(node_id), tolerance, check_delay);
                % perfect delay correction 
                % constant_delay_1_hop(node_id,k) = delay(node_id);
            end

            % query data from database (ie, propagate stored data to 
            % current V time to use for CCF)
            [N, ~] = DQF(Database_Stored{node_id}, V(node_id, :), "update", delay_confirmed(node_id));%, H(node_id,:), rate(node_id,:), dH(node_id, :), CC, delay_confirmed(node_id));
            
            %Database_Stored{node_id}(:,6) = Database_Stored{node_id}(:,6)+1;
            
            % filter received data
            [N_filtered, outlier_check(node_id)] = CCF_filter(filter_type, N, 0, nb_occurences_per_ID);

            % anotate database - data not filtered out 
            Database_Stored{node_id} = database_anotate(Database_Stored{node_id}, N, V, N_filtered);
            
            % percentage data NOT-filtered
            percentage_filtered(node_id,k) = length(N_filtered)/length(N);
            
            % update rate                   
            [rate(node_id, k), C(node_id), D(node_id), C_rate(node_id)] = CCF_update(delta_H(node_id), N_filtered, V(node_id,k), rate(node_id, k), C(node_id), D(node_id), C_rate(node_id), alpha, overlap(node_id), database_check, outlier_check, eta, gamma, smallest_rate);                      
            
            % change CC cycle for this node if update does not correspond
            % to initial CC
            if change_update_pattern(node_id)
                CC_adapted(node_id) = delta_H(node_id);
            end
            % calculate next update trigger
            H_update(node_id) = H_trigger(H(node_id, k), CC_adapted(node_id), dH(node_id), "future");


            last_update(node_id)=H(node_id,k);
            updates(node_id,k)=1;
        end
        
        % delays 
        stored_delay_confirmed(:,k)=delay_confirmed;
        stored_delay_rate(:,k)=delay_rate;

        % Convergence analysis/Calculate normalised node difference
        % only check non-disruptive nodes
        % error in time
        relative_time_diff_k = node_diff(V(i_normal,k), max_init_time, min_init_time);
        relative_time_diff = cat(2,relative_time_diff, relative_time_diff_k);
        
        
        % error in rate
        relative_rate_diff_k = node_diff(absolute_rate(i_normal,k), max_init_rate, min_init_rate) ;
        relative_rate_diff = cat(2,relative_rate_diff, relative_rate_diff_k);

        % simulate while max time step not exceeded 
        if k>max_time
            % maximum errors between non-disruptive nodes 
            % in time & rate
            max_error_time = max_difference(V, i_normal);
            max_error_rate = max_difference(absolute_rate, i_normal);
            break
        end
    end  
    
    
    % when debugging/correcting
    % opinion_plots_single(V, t,0, ...
    %                 0, i_dis, [true false false], ...
    %                 10^(-9),"mono", "Time")
    % 
    % opinion_plots_single(absolute_rate, t, 0, ...
    %                 0, i_dis, [true false false], ...
    %                 10^(-9),"mono", "Rate")
 
    % figure()
    % plot(repmat(t, [n, 1]).',V.','-')
    % figure()
    % plot(repmat(t, [n, 1]).',absolute_rate.','-')
   

end

