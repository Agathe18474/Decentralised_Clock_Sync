function [data, ID] = DQF(Database, V, purpose, delay)
%DQF(Database, V, purpose, delay, H,  rate, dH, CC)
    % Database Querying Function - extract data of use, propagate it to 
    % current V-time
    % Inputs: 
    %   Database                database of node's data
    %   V                       node's current V-time
    %   purpose                 purpose of DQF
    %                               > "broadcast"
    %                               > "update"
    % Outputs: 
    %   data

    % CHANGED - ADDED DELAYS!!!!

    N =  Database(:,1); 
    ID =  Database(:,2);
    hop_check = Database(:,5);
    filtered_out_check = Database(:,6);
    %H_at_reception = Database(:,3);
    
    %CC_usage = Database(:,6);
    switch purpose 
        case "broadcast"
            % only broadcast 1-hop neighbour data
            %idx_data_to_broadcast = Database.hop_check ==true;    
            idx_data_to_broadcast = hop_check ==1;  
            
            % flagged data (anotated) = not broadcasted
            % idx_1_hop = hop_check ==1;  
            % idx_filtered_check = filtered_out_check == 1;
            % idx_data_to_broadcast = and(idx_1_hop,~idx_filtered_check);


            
            if ~isempty(idx_data_to_broadcast)
                data_output = N(idx_data_to_broadcast);
                ID = ID(idx_data_to_broadcast);
                V_at_reception = Database(idx_data_to_broadcast,4);
                
            else
                V_at_reception = [];
                data_output = [];
                ID = [];
            end
        case "update"
            % all data received 
            %data_output = N;
            V_at_reception = Database(:,4);
            H_at_reception = Database(:,3);
            data_output = N;
            ID = [];
        otherwise
            error('DQF Purpose not recognised')
    end

    if ~isempty(data_output)
        % propagation time is diff between arrival time & current
        % V-time
        switch purpose  
            case "broadcast"
                data = data_output-(V_at_reception+delay);
                data = data+V(end);
            case "update" 
                data_output = N-(V_at_reception+delay);
                data = [data_output,H_at_reception];
                       
        end
    else 
        if purpose == "broadcast"
            data = data_output;
        else
            data = double.empty([0 2]);
        end
        
    end
 
    % propagate data forward
    %data = propagator(data, propag_time);
end