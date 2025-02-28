function [Database_Filtered, nb_occurences] = database_filter(Database_Received, self_id, step, sampling, nb_occurences_2_hop_Step_1)
%function [Database_Filtered, nb_occurence_ID, nb_occurence_ID_2_hop] = database_filter(Database_Received, self_id, step, sampling, nb_occurences_2_hop_Step_1)
    
% Filter database of received values
    % Keep 1-hop transmission over 2-hop transmission from same node
    % Keep random or median 2-hop transmission received from node 
    % Inputs: 
    %   Database_Received           database of received values
    %   self_id                     id of node database belongs to
    %   step                        step of filtering 
    %   sampling                    how is 2-hop data selected? 
    %                                   > "random": pick 1 data pt at rdm
    %                                       for each ID
    %                                   > "median": pick median data pt of
    %                                       each ID
    %   nb_occurences_2_hop_Step_1  nb data pts for each 2-hop ID from prev
    %                                   filtering step
    % Ouputs: 
    %   Database_Filtered           filtered database 
    %   nb_occurence_ID             nb data pts received for each ID
    
    % extract values from table for ease of use
    ID = Database_Received(:,2);
    
    % remove data with same self id
    Database_Received(ID == self_id, :)=[];

    % extract IDs and hop_check from table for ease of use
    ID = Database_Received(:,2);
    hop_check = Database_Received(:,5);
    
    % keep 1-hop data when available over 2-hop data 
    % keep only latest 1-hop data
    data_1_hop_ID = ID(hop_check==1);
    Database_1_hop = Database_Received(hop_check==1, :);
    [~, idx_flip, ~] = unique(flip(data_1_hop_ID));
    idx_flip_2 = length(data_1_hop_ID)-idx_flip+1;
    Database_1_hop = Database_1_hop(idx_flip_2, :);

    nb_occurence_ID_1_hop = ones(size(Database_1_hop,1),1);

    %[~, sorted_id] = sort(Database_Received(hop_check==1, :);
    
    % find redundant rows - ie rows of 2-hop data for nodes where 1-hop
    % data is available
    data_2_hop_ID =  ID(hop_check==0);
    redundant_ID = data_1_hop_ID(ismember(data_1_hop_ID, data_2_hop_ID));
    redundant_rows = and(ismember(ID, redundant_ID),hop_check == 0);
    
    % remove redundant rows from database
    Database_Received(redundant_rows,:) = [];
    
    % re-extract values from table for ease of use
    ID = Database_Received(:,2);
    hop_check = Database_Received(:,5);
    data_2_hop_ID = ID(hop_check == 0);
     
    % find node IDs for which multiple 2-hop data has been received
    U = unique(data_2_hop_ID);
    
    % if no 2-hop data
    if isempty(U) 
        Database_2_hop_filtered = Database_Received(hop_check == 0, :);
        nb_occurence_ID_2_hop=[];
        %nb_occurence_ID = nb_occurence_ID_1_hop;
    else
        Database_2_hop = Database_Received(hop_check == 0, :);
        % select data to keep 
        % ensures that node values are not selected from a single transfer
        % ie, that they have not all been transfer through the same
        % neighbour

        % if no duplicates (no data from same ID)
        if length(unique(data_2_hop_ID))==length(data_2_hop_ID)
            Database_2_hop_filtered = Database_2_hop;
            nb_occurence_ID_2_hop = ones(length(data_2_hop_ID),1);
        else 
            if step == "Step 1"
                % sort 2-hop data by ID
                [sorted_ID, idx_sorted] = sort(data_2_hop_ID);
                Database_2_hop_sorted = Database_2_hop(idx_sorted, :);

                % create indexes to sample data from sorted data array
                % check when transition to new ID
                % eg: ID = [1, 1, 2, 2, 3, 3, 3], idx_diff_ID = [2, 4]
                idx_diff_ID = find(diff(sorted_ID));
                % create start & end of ID samples 
                % eg: start: [1, 3, 5]
                % eg: end: [2, 4, 7]
                start_idx_sample = [1; idx_diff_ID+1];
                end_idx_sample = [idx_diff_ID;length(sorted_ID)];
                
                % unique IDs
                [edges, ~,~] = unique(sorted_ID);   
                % # time IDs appear in Database                 
                nb_occurence_ID_2_hop = histcounts(data_2_hop_ID(:), [edges;inf]).';
                % nb_occurence_ID_2_hop(:,2) = edges;
                % nb_occurence_ID = [nb_occurence_ID_1_hop;nb_occurence_ID_2_hop(:,1)];

                if sampling == "median"
                    % assume that values are already sorted
                    % which they are as they are placed in order of arrival
                    % find mid point of start & end ID sample
                    % ie, select median of each ID sample since values
                    % sorted
                    mid_pt = floor((start_idx_sample+end_idx_sample)./2);

                    % val_sorted_2 = accumarray(sorted_ID, val_sorted, [], @sort);                    
                    % [median_idx2,~] = ismember(vals_sorted,val_sorted_2(mid_pt));    
                    % Database_2_hop_filtered2 = Database_2_hop_sorted(median_idx2, :);

                    % select mid-point of each ID data
                    Database_2_hop_filtered = Database_2_hop_sorted(mid_pt, :);

                elseif sampling == "random"
                    % randomly sample idx for each ID
                    % ie, randomly select 1 datapoint for each ID 
                    sampled_idx = start_idx_sample + round((end_idx_sample-start_idx_sample).*rand(length(start_idx_sample), 1));   
    
                    Database_2_hop_filtered = Database_2_hop_sorted(sampled_idx, :);
                end
            elseif step == "Step 2"
                % select latest data
                H_at_reception_2_hop = Database_2_hop(:,3);
                [~, idx_sorted] = sort(H_at_reception_2_hop, "descend");
                Database_2_hop_Sorted_Reception = Database_2_hop(idx_sorted,:);
                [~,idx_unique] = unique(Database_2_hop_Sorted_Reception(:,2));
                Database_2_hop_filtered = Database_2_hop_Sorted_Reception(idx_unique, :);
                
                % match IDs which were from Database_Received
                match_ID_Step_1 = ismember(Database_2_hop_filtered(:,2), nb_occurences_2_hop_Step_1(:,2));
                
                % some 2-hop data from given ID removed at benefit of 1-hop 
                % transmission for same ID 
                remove_new_1_hop_data = ismember(nb_occurences_2_hop_Step_1(:,2),Database_2_hop_filtered(:,2));
                
                % nb_occurences for IDs which were from received kept
                %nb_occurence_ID_2_hop = ones(size(Database_2_hop_filtered, 1),1);
                %nb_occurence_ID_2_hop(match_ID_Step_1) = nb_occurences_2_hop_Step_1(remove_new_1_hop_data,1);

            end
                  

    
        end
    end        
    % output cleaned database
    Database_Filtered = [Database_1_hop; Database_2_hop_filtered];
    nb_occurences = [nb_occurence_ID_1_hop; nb_occurence_ID_2_hop];

    
    % TEST
    % keep odd number of data!
    % nb_data_points = size(Database_Filtered, 1);
    % if rem(nb_data_points,2)==0
    %     random_nb = randsample(nb_data_points, 1);
    %     Database_Filtered(random_nb,:)=[];
    % end
    % elseif step == "Step 2"
    %     if length(unique(Database_Received(:,2)))~=length(Database_Received(:,2))
    %         H_at_reception = Database_Received(:,3);
    %         [~, idx_sorted] = sort(H_at_reception, "descend");
    %         Database_Sorted_Reception = Database_Received(idx_sorted,:);
    %         [~,idx_unique] = unique(Database_Sorted_Reception(:,2));
    %         Database_Filtered = Database_Sorted_Reception(idx_unique, :);
    %     else
    %         Database_Filtered = Database_Received;
    %     end
    %end
end