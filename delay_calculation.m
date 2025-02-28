function [delay, delay_rate, check_delay] = delay_calculation(Database, dH, delay_past, alpha, delay_rate_past, min_nb_sources)
    % Estimate tof delays 
    % Inputs:
    %   Database                databse of received data
    %   dH                      elapsed H-time since last estimate
    %   delay_past              past delay
    %   alpha                   maximum smoothing factor
    %   delay_rate_past         past delay rate
    %   min_nb_sources          min # 1-hop connections
    % Outputs:
    %   delay                   new delay
    %   delay_rate              new delay rate
    %   check_delay             enough data to implement delay correction

    % received data from 1-hop neighbours
    idx_1_hop = Database(:,5)==1;

    % errror
    error_1_hop = Database(idx_1_hop, 1)-Database(idx_1_hop, 4);
    delay =  mean(error_1_hop);
    
    % delay rate
    delay_rate_new = (delay-delay_past)/dH;
    
    % smooth delay rate
    delay_rate = single_exponential_smoothing(alpha, delay_rate_new, delay_rate_past);

    check_delay = false;
    if sum(idx_1_hop)>min_nb_sources
        check_delay = true;
    end
end