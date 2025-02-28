function delay_confirmed = delay_confirmation(delay_confirmed, delay_new, delay_rate, tolerance, check_enough_data)
    % Confirm delay if meet conditions:
    %       > delay rate is stable, within tolerance
    %       > enough data points used to calculate delay
    %       > delay is below max range
    % Inputs: 
    %   delay_confirmed     previous confirmed delay
    %   delay_new           new delay
    %   delay_rate          new delay rate
    %   tolerance           tolerance for rate
    %   check_enough_data   logical check enough datapoints
    % Outputs: 
    %   delay_confirmed     confirmed delay


    % max ISL range in km
    max_range = 20000;
    % speed of light in km/s
    c = 299792.458;
    max_delay = max_range/c;
    % check if delay is negative & below max_delay
    check_within_range = max_delay>abs(delay_new) && delay_new<0;

    % check if rate is within tolerance
    check_rate_tolerance = abs(delay_rate)<=tolerance;
 
    if check_rate_tolerance && check_enough_data && check_within_range
        delay_confirmed = 0.8*delay_new;
    end
        
end