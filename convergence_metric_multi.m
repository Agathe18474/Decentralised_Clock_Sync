function [convergence_metric,convergence_metric_min, convergence_metric_max]  = convergence_metric_multi(convergence_metric_single)
    % Calculate convergence metric for a SINGLE run 
    % Inputs: 
    %   relative_diff               relative difference between nodes
    %   zero_val                    relative zero val acting as floor
    % Outputs: 
    %   convergence_metric          mean of monte carlo convergence metrics
    %   convergence_metric_min      min of monte carlo convergence metrics
    %   convergence_metric_max      max of monte carlo convergence metrics

    convergence_metric = mean(convergence_metric_single, 1, "omitnan");
 
    convergence_metric_min = min(convergence_metric_single, [], 1) ;
    convergence_metric_min(1, convergence_metric_min==0) = convergence_metric(1, convergence_metric_min==0);
    convergence_metric_max = max(convergence_metric_single, [], 1) ;

end