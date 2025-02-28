function [N_filtered, outlier_check] = CCF_filter(filter_type, N, V, nb_occurences_per_ID)
    % Concensus Calculation Function - filtering of incoming data 
    % Remove extreme values 
    % Inputs: 
    %   filter_type                 type of filter used to remove outliers
    %                                   > ODDI-C
    %                                   > MSR
    %                                   > Hampel
    %                                   > Dixon
    %                                   > Mean-based
    %   N                           neighbour values (1- & 2-hops)
    %   V                           node's current V-time
    %   nb_occurences_per_ID        nb data pts received from that ID
    % Outputs: 
    %   N_filtered                  filtered neighbour values (1- & 2-hops)
    %   outlier_check               check (true/false) if node's own V-time
    %                               is an outlier 
    

     
    % Basic 
    % N = N(:,1);
    
    % TEST 
    % Weight calculated based on nb data pts received from that ID 
    weights = 0.8.^(nb_occurences_per_ID-1);
    N = weights.*N(:,1); 
    

    if isempty(N)
        N_filtered = N;
        outlier_check = false;
    else
        switch filter_type
            case "ODDI-C"
                [N_filtered, outlier_check] = filter_ODDI_C(N, V);                
            case "MSR"
                N_filtered = filter_MSR(N, i_dis, V);
            case "Hampel"
                filter_val = 3;
                [N_filtered, outlier_check] = filter_Hampel(N, V, filter_val);
            case "Dixon"
                [N_filtered, outlier_check] = filter_Dixon(N, V);
            case "Mean-based"
                N_filtered = N;
                outlier_check = false;
        end
    end
end