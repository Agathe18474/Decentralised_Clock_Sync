function dis_param = disruptor_parameters(func, max_init, min_init, coordinated, max_time_step)
    % byzantine node parameters for given type of func chosen
    % Inputs: 
    %   func            type of node behaviour  
    %                   sine_large_amp: large amplitude sine wave
    %                   sine_small_amp: small amplitude sine wave
    %                   const: constant 
    %                   noise: random 
    %   coordinated     whether nodes show coordinated behaviour (boolean)
    % Output:
    %   byz_param       byzantine node parameter array 
    %                   [period, shift_x,   shift_y]xNb nodes
    
    % find sine type disruptive nodes 
    sine_dis_nodes = contains(func, "sine"); 
    const_dis_nodes = contains(func, "const"); 
    % initialise byzantine parameter array 
    dis_param = zeros(length(func), 4);

    if sum(sine_dis_nodes)+sum(const_dis_nodes)==0
        return
    end
    
    % i nodes coordinated -  all sine type nodes have the same parameters 
    if (nnz(sine_dis_nodes)>0) && coordinated 
        angular_freq = 4*rand;
        amplitude = 0.8*rand;
        shift_x = 5*rand;
        shift_y = max_init - (max_init-min_init)*rand;
        dis_param(sine_dis_nodes, :) = [angular_freq, amplitude, shift_x, shift_y];
        
    elseif (nnz(sine_dis_nodes)>0) && coordinated == false
        param = [];

        % loop through all sine type nodes
        for i=1:nnz(sine_dis_nodes)
            angular_freq = 4*rand;
            shift_x = 5*rand;
            shift_y = rand*((max_init - (max_init-min_init)*rand)-(max_init+min_init)/2);
            amplitude = 0.8*rand;
            param_new = [angular_freq, amplitude, shift_x, shift_y];
            param = cat(1,param,param_new);
        end

        dis_param(sine_dis_nodes, :) = param;
    end
    if nnz(const_dis_nodes)>0
        for i=1:nnz(const_dis_nodes)
            m(i,1) = 2*(max_init-min_init)*(-1+2*rand)/(max_time_step);   
        end
         dis_param(const_dis_nodes, 1)=m;
    end
end