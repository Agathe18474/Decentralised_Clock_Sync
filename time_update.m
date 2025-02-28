function val_future = time_update(val_current, rate, step) 
    % predict update based on exponential average 
    % Inputs:
    %   val_current             current value
    %   step                    step size between current and past val
    %   rate                    rate of update
    %   alpha                   exponential factor
    % Outputs: 
    %   val_future              future value

    % predicted value 
    val_future = rate*step + val_current;
end