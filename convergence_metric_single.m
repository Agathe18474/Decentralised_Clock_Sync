function convergence_metric_single = convergence_metric_single(relative_diff, zero_val)
    % Calculate convergence metric for a SINGLE run 
    % Inputs: 
    %   relative_diff               relative difference between nodes
    %   zero_val                    relative zero val acting as floor
    % Outputs: 
    %   convergence_metric_single   convergence metric for single run

    convergence_metric_single = relative_diff(1,:)./relative_diff(1,1);
    convergence_metric_single(convergence_metric_single<zero_val)=zero_val;   

end