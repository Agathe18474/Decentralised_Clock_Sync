function dis_node_val = disruptor_node_val(func, param, max, min, t)
    % byzantine node value at given type step based on type of func chosen
    % Inputs: 
    %   func            type of node behaviour  
    %                   sine_large_amp: large amplitude sine wave
    %                   sine_small_amp: small amplitude sine wave
    %                   sine_very_small_amp: very small amplitude sine wave
    %                   const: constant 
    %                   noise: random 
    %   param           parameters for sine functions 
    %                   [period, shift_x,   shift_y]
    %   max             maximum possible value that cannot be exceeded
    %   min             minimum possible value that cannot be exceeded
    %   t               current time
    % Output:
    %   dis_nod_val     value of byzantine node at time t

    % values 
    avg = (max+min)/2;
    range = max-avg;
    
    switch func
        case "sine_wave"
            dis_node_val = (avg+param(1,4)) + (param(1,2)*range)*sin( (param(1,1)*t)+param(1,3) ) + param(1,1)*t;
            %dis_node_val = (avg+param(1,3)) + (0.8*range)*sin( (param(1,1)*t)+param(1,2) );
        case "sine_small_amp"
            dis_node_val = (avg+param(1,3)) + (0.3*range)*sin( (param(1,1)*t) + param(1,2) );
        case "sine_very_small_amp"
            dis_node_val = (avg+param(1,3)) + (0.1*range)*sin( (param(1,1)*t) + param(1,2) );
        %case "const"
         %  dis_node_val = (avg + range*0.3 + 0.3)*ones(1,t(end));
        case "lin_const"
            %dis_node_val = avg + rand*0.5*t;
            dis_node_val = avg + param(1,1)*t;
        case "noise"
            a = rand;
            if a>0.5
                dis_node_val = avg + rand(1, t(end))*0.9*range;
            else 
                dis_node_val = avg - rand(1, t(end))*0.9*range;
            end
    end

end