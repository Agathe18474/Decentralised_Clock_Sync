function events = event_triggering(V, V_last_broadcast)
    % Event triggering function 
    % Inputs: 
    %   V                   current V-time
    %   V_last_broadcast    V-time at instance of last broadcast
    % Outputs: 
    %   events              boolean of nodes having exceeded threshold
    
    thresh_1 = 0.05;
    thresh_2 = 0.2;

    % difference between latest broadcast and current value
    %error_broadcast = abs(V-V_last_broadcast);

    error_broadcast = abs(V-V_last_broadcast)./V;

    events = error_broadcast.^2>thresh_1;
    
    % WORKS FOR RINGS!!!
    %events = error_broadcast>thresh_1*exp(-thresh_2*V_last_broadcast);
    %events = error_broadcast>thresh_1*exp(-thresh_2*error_broadcast);


end