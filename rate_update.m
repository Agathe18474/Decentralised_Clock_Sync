function smoothed_rate_current = rate_update(val_current, val_past, step_size, smoothed_rate_past, alpha)
    % predict update based on exponential average 
    % Inputs:
    %   val_current             current value
    %   val_past                past value
    %   smoothed_rate_past      past weighted average of rate
    %   alpha                   exponential factor
    % Outputs: 
    %   smoothed_val_current    new weighted average
    % Function called: 
    %   single_exponential_smoothing

    % current rate of growth for this step
    rate_current = (val_current-val_past)/step_size;
    
    % single exponential smoothing
    smoothed_rate_current = single_exponential_smoothing(alpha, rate_current, smoothed_rate_past);

end