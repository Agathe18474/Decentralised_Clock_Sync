function [rate, C, D, C_rate] = CCF_update(delta_H, N_filtered, V, rate_past, C_past, D_past, C_rate, alpha, overlap, database_check, outlier_check, eta, gamma, smallest_rate)
    % Consensus Cycle Function - Update virtual rate
    % Inputs:   
    %   dH_CC                   consensus cycle in terms of dH
    %   N_filtered              filtered neighbour errors (1- & 2-hops)
    %   V                       current V-time 
    %   rate                    current V-rate
    %   C_past                  consensus time at prev CC
    %   D                       Past average of filtered errors
    %   C_rate                  consensus rate at prev CC
    %   alpha                   maximum smoothing factor 
    %   overlap                 amount of overlap between neighbour's IDs
    %                              between CCs
    %   database_check          check if database is empty [true/false]
    %                               if not - implies node is disconnected
    %                               since the last PC (Prune Cycle)
    %   outlier_check           check if node's V is outlier [true/false]
    %   eta                     learning factor of error in time
    %   gamma                   learning factor of error in rate
    %   smallest_rate           min value new rate allowed to take on
    % Outputs: 
    %   rate                    V-rate of current CC
    %   C                       consensus time of current CC
    %   D                       average of filtered errors of current CC
    %   C_rate                  C-rate of current CC

    % Functions used: 
    %   rate_update
    %   single_exponential_smoothing


    % damper factor for C-rate
    %damper = (alpha*overlap)^(CC);
    damper = (alpha*overlap);
    %damper = alpha;
    %damper = overlap; 
        
    if  C_past==0 && database_check
        % if database is empty AND no data from prev CC
        eta = 0;
        gamma = 0;  
        
        D = 0;
        C = 0;

        % no data available - no error 
        error_t = 0;
        error_rate= 0;

    elseif C_past==0 && ~isempty(N_filtered)
        % if database NOT empty BUT no C-time from prev CC

        % calculate C-time
        % if node's own V-time is outlier (outlier check) - ignore it for 
        % the C-time calculation
        if outlier_check
            %C = V+mean(N_filtered);
            D = mean(N_filtered);
        else
            %C = V+mean([N_filtered; 0]);
            D = mean([N_filtered; 0]);
        end

        D = single_exponential_smoothing(damper, D, D_past);
        C = V+D;

        % errors
        % error in time only  - error in rate unavailable 
        error_t = (C-V)/delta_H;
        
        error_rate=0;

    elseif  isempty(N_filtered) || database_check
        % if all values from database have been filtered out, ie node's
        % V-time is ON the median.

        % C is the node's V-time
        D = 0;
        C = V;

        % calculate C-rate 
        % C-rate based on current C-rate and C-rate of prev CC
        % prevents drastic changes in rate when node reconnects or receives
        % new data
        % smoothing factor dictates if preference given to past or new data
        C_rate_past = C_rate;

        %C_rate = rate_update(C, C_past, delta_H,  C_rate_past, damper);
        C_rate = rate_update(C, C_past, delta_H,  C_rate_past, 1);
        
        % errors
        error_rate = (C_rate - rate_past);
        error_t = 0;
       
    else 
        % if database not empty & data from prev CC

        % calculate C-time
        % if node's own V-time is outlier - ignore it for the C-time
        if outlier_check
            %C = V+mean(N_filtered);
            D = mean(N_filtered);
        else
            %C = V+mean([N_filtered; 0]);
            D = mean([N_filtered; 0]);
        end

        D = single_exponential_smoothing(damper, D, D_past);
        C = V+D;

        % calculate C-rate 
        % C-rate based on current C-rate and C-rate of prev CC
        % prevents drastic changes in rate when node reconnects or receives
        % new data
        % smoothing factor dictates if preference given to past or new data
        C_rate_past = C_rate;

        %C_rate = rate_update(C, C_past, delta_H,  C_rate, damper);
        C_rate = rate_update(C, C_past, delta_H,  C_rate_past, 1);
          
        % errors
        error_t = (C-V)/delta_H;
        error_rate = (C_rate - rate_past);        
    end
    
    % total error in rate
    if sign(error_t) == sign(error_rate) 
        % if error in time & rate have the same sign
        % take biggest error 
        error =  alpha*sign(error_t)*max(eta*abs(error_t), gamma*abs(error_rate));
    else
        % otherwise combine the two together 
        error = alpha*(eta*error_t + gamma*error_rate);
    end 

    % update rate 
    %rate = rate + error;
    rate = single_exponential_smoothing(damper, rate_past+error, rate_past);   
               
    % minimum rate
    if rate<smallest_rate
        rate = smallest_rate;
    end

    
end