function [conn_val_filt, outlier_check] = filter_Hampel(conn_val, node_val_id, filter_val)
    % Hampel outlier filter
    %
    % Inputs:
    %   conn_val        connected nodes' values 
    % Outputs:
    %   conn_val_filt   connected nodes' filtered out values
    %   outlier_check   check if node is an outlier   

    % median
    med = median(conn_val);
    % median absolute deviation
    MAD = median(abs(conn_val-med));    
    % scaled MAD
    MAD_scaled = MAD/0.6745 ;
    % connected values z-score (median based)
    conn_z = abs((conn_val-med)./MAD_scaled);

    z_id = abs((node_val_id-med)./MAD_scaled);
    outlier_check = false;
    if z_id>3
        outlier_check = true; 
    end
    
    % Hampel identifier, values outside of +/- 3*NMAD 
    conn_val_filt = conn_val(conn_z<=filter_val);
    
end