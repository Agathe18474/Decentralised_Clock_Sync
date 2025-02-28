function Database_Appendage = database_appender(data_broadcast, ID_broadcast, V, H, tofs, V_broadcastor, ID_broadcastor, outlier_check)
    % Formatting of received data into database 
    % Inputs: 
    %   data_broadcast          data from broadcastor node
    %   ID_broadcast            ID of values forwarded by broadcastor
    %                               node
    %   V                       V-time of receptor node
    %   H                       H-time of receptor node
    %   tof                     transmission delay between broadcasting &
    %                           receptor nodes based on tof (time of flight)
    %   V_broadcastor           V-time of broadcastor node
    %   ID_broadcastor          ID of broadcastor node
    % Ouputs: 
    %   Database_New            new data formatted for database  

    outlier_check = false;
    
    % length of appended data
    if outlier_check==false
        new_data_length = length(data_broadcast)+1;
    
        % 2-hop neighbours' V-time & broadcastor V-time
        N = [data_broadcast; V_broadcastor];  
        % ID of N-times
        ID = [ID_broadcast; ID_broadcastor];
    
        % H- & V-time at arrival of data    
        H_at_reception = (H + tofs(1))*ones(new_data_length,1);
        V_at_reception = (V + tofs(2))*ones(new_data_length,1);
    
        % 1-hop check == false as this is 2-hop data
        hop_check = zeros(new_data_length,1);    
        hop_check(end) = 1;

    else
        new_data_length = length(data_broadcast);

        % 2-hop neighbours' V-time & broadcastor V-time
        N = data_broadcast;  
        % ID of N-times
        ID = ID_broadcast;
    
        % H- & V-time at arrival of data    
        H_at_reception = (H + tofs(1))*ones(new_data_length,1);
        V_at_reception = (V + tofs(2))*ones(new_data_length,1);
    
        % 1-hop check == false as this is 2-hop data
        hop_check = zeros(new_data_length,1);   
    end

    % CC usage - for how many CC has the data been used
    CC_usage = zeros(new_data_length,1);
    
    % New data for database  
    Database_Appendage = [N, ID, H_at_reception, V_at_reception, hop_check, CC_usage];
end