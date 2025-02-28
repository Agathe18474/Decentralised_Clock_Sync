function H_trigger = H_trigger(H_current, cycle_length, dH, type)
    % Create trigger time past which subsequent function is activated
    % Inputs:
    %   H_current           current H-time
    %   cycle_length        length of cycle after which trigger activates
    %   type                type of trigger 
    %                           "future" - trigger is in the future
    %                           "past" - trigger is in the past (to check
    %                               when data no longer relevant)
    % Outputs:
    %   H_trigger           trigger H-time
    
    % trigger time
    if type == "future"
        %H_trigger = H_current + cycle_length;
        H_trigger = H_current + 0.999*cycle_length;
    elseif type == "past"
        H_trigger = H_current - cycle_length;
        %H_trigger = H_current - 1.011*cycle_length;
    end
end