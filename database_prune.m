function Database_Stored_pruned = database_prune(Database_Stored, H_prune)
    % Prune database of older, no longer relevant values
    % Data is pruned when collected before the start of the Prune Cycle
    % (PC) indicated by the H-time H_prune
    % Inputs: 
    %   Database_Stored             database of current stored values
    %   H_prune                     H-time before which data is pruned
    % Ouputs: 
    %   Database_Stored_pruned      database of current stored values
    %                                   pruned of older data
    
    % 
    H_at_reception = Database_Stored(:,3);
    

    % prune older data
    Database_Stored_pruned = Database_Stored;
    prune_idx = H_at_reception < H_prune;
    Database_Stored_pruned(prune_idx, :) = [];

    %prune_idx = Database_Stored_pruned.H_at_reception < H_prune;
    %Database_Stored_pruned(prune_idx, :) = []; 

end