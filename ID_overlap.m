function ID_overlap = ID_overlap(ID_past, ID_current)
    % Overlap calculation - extent to which current dataset has grown or
    % shrinked from previous timeseto
    % Inputs: 
    %   ID_past                 ID of neighbours' nodes (1- & 2-hops) from  
    %                               previous CC
    %   ID_current              ID of neighbours' nodes (1- & 2-hops) from  
    %                               current CC
    % Outputs: 
    %   ID_overlap              percentage overlap between ID past & new
    
    % have nodes been lost from previous CC? what is the % of nodes from
    % prev CC which are still connected?
    overlap_backwards = length(intersect(ID_past, ID_current))/length(ID_past);

    % have nodes been added from previous CC? what is the % of nodes from
    % current CC which are new connections?
    overlap_forwards = length(intersect(ID_past, ID_current))/length(ID_current);
    
    % select smallest
    ID_overlap = min(overlap_backwards, overlap_forwards);
end