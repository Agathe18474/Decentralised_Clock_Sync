function Database_Anotated = database_anotate(Database_Stored, N, V, N_filtered)
    % Filter database of received values
    % Keep 1-hop transmission over 2-hop transmission from same node
    % Keep last 2-hop transmission received from node (based on V-time at
    % reception)
    % Inputs: 
    %   Database_Received           database of received values
    % Ouputs: 
    %   Database_Filtered           filtered database

    %N_filtered_Hampel = CCF_filter("Hampel", N, V);

    Database_Anotated = Database_Stored;
    %ID_anotate = ismember(N(:,1),N_filtered_Hampel);
    ID_anotate = ismember(N(:,1),N_filtered);

    %Database_Anotated(~ID_anotate,6) = Database_Stored(~ID_anotate,6)+1;

    % if anotate ~ID_anotate, flag data filtered OUT
    % if anotate ID_anotate, flag data that has gone through filtering

    % data that is flagged - not broadcasted 
    Database_Anotated(:,6) = ~ID_anotate;


end