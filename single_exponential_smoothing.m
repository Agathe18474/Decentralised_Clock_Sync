function smoothed_val_current = single_exponential_smoothing(alpha, val_current,smoothed_val_past)
    % Exponential weighted average
    % Inputs:
    %   alpha                   exponential factor
    %   val_current             current value
    %   smoothed_val_past       past weighted average
    % Outputs: 
    %   smoothed_val_current    new weighted average 

    smoothed_val_current = alpha*val_current + (1-alpha)*smoothed_val_past;

end