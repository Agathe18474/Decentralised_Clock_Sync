function [updating_nodes, forced_update] = update_trigger(H, H_update, Database_size, Database_nb_sources, min_size, max_size, min_nb_sources)
    % check which nodes meet updating conditions 
    % size of received database must be: 
    %       > larger than min_size & at least 3 1-hop connections
    %       > smaller than max_size & delta_H<CC 
    % Inputs:
    %   H                       current H-times
    %   H_update                H-times at which update planned
    %   Database_size           size received database, ie # entries
    %   Database_nb_sources     nb of 1-hop connections
    %   min_size                min # entries to allow update
    %   max_size                # entries after which force update
    %   min_nb_sources          min # 1-hip connections to allow update
    %                               when data no longer relevant)
    % Outputs:
    %   updating_nodes          logical array of updating nodes
    %   forced_update           logical array of forced updates
    
    % unless meet min condition, do NOT update
    % larger than min_size & at least 3 1-hop connections 
    min_size_check = and(Database_size>min_size, Database_nb_sources>min_nb_sources);

    % smaller than max_size & delta_H<CC
    %max_size_check = and(Database_Received_size>max_size_database, delta_H<CC);
    %max_size_check = and(Database_size>max_size, Database_nb_sources>min_nb_sources);
    max_size_check = and(Database_size>min_size, Database_nb_sources>10);

    %database_size_check = and(Database_Received_size>max_size_database, delta_H<CC);
    trigger = H>=H_update;
    min_conditions = and(trigger, min_size_check); 
    

    % update if meet min conditions OR if max conditions have been exceeded
    % and no update was set to occur
    %updating_nodes = or(H(:,k)>=H_update,database_size_check);
    updating_nodes = or(min_conditions,max_size_check);
    
    forced_update = false(length(H),1);
    if any(and(~min_conditions,max_size_check))
        forced_update(max_size_check) = true;
    end
end