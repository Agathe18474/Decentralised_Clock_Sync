function [A_time_step, tof_time_step] = time_step_adajacency(A, A_type, k, tot_sat, tof)
    % DESCRIPTION:
    %   Extract adjancency matrix of current time step. 
    %
    % INPUTS:
    % A             Adjacency variable
    % A_type        Class of Adjacency variable
    % k             Current time step 
    % OUTPUTS: 
    % A_time_step   Adjancency matrix of current time step

    if contains(A_type, 'logical')
        % satellite contact plans
        % extract each time steps plan 
        Counter_temp = k-floor((k-1)/size(A,2))*size(A,2);
        A_time_step = reshape(A(:,Counter_temp), tot_sat, tot_sat);
        if sum(sum(tof))==0
            tof_time_step = zeros(tot_sat);
        else
            tof_time_step = reshape(tof(:,Counter_temp), tot_sat, tot_sat);
        end
    elseif contains(A_type, 'cell')
        % satellite network contained in structure
        % extract current time step's matrix 
        % when reach end of structure, cycle back to start
        Counter_temp = k-floor((k-1)/length(A))*length(A);
        A_time_step = cell2mat(A(Counter_temp));
        A_time_step = full(A_time_step);
        A_time_step = A_time_step+A_time_step.';

        % create random distances 
        % based on ISL = 4000km 
        % create rdm distribution between 1000 & 4000
        % speed of light      
        %c = 299792.458;
        %tof_random = 3000*rand(tot_sat) + 1000;
        %tof_time_step = tof_random./c;
        %tof_time_step = tof_time_step.*A_time_step;
        %tof_time_step = 0.0077*(A_time_step~=0);
        
        % just take in fixed input data
        tof_time_step = tof;
    elseif contains(A_type, 'double') && size(A,3)>1 
        % rolling matrix
        % extract current time step's matrix 
        % when reach end of structure, cycle back to start
        Counter_temp =k-floor((k-1)/size(A,3))*size(A,3);
        A_time_step = full(A(:,:,Counter_temp));
    elseif contains(A_type, 'double') 
        % static network 
        A_time_step = A;
        tof_time_step =  tof;
    else
        error('Invalid adjacency matrix type')
    end

end