function [conn_val_filt, outlier_check] = filter_ODDI_C(conn_val, node_id_val)
% Opinion Dynamics-inspired Disruption-tolerant Consensus (ODDI-C) Filter
%
% DESCRIPTION: 
%   Filters values out if their median z-score are smaller than the node's 
%   own median z-score. Ensures that the more of an outlier a node is, 
%   the more willing it will be to account for other opinions while a 
%   node fitting in with other nodes will stay put.
%
%   Inspired by Deffuant-Weisbuch opinion model & adapted using 
%   Robust Statiscal Approaches (median, MAD, and scaled MAD)
%
% REFERENCES:
%   [1] P. Sobkowicz, "Extremism without extremists: Deffuant model 
%   with emotions", Front Phys. Interdisciplinary Phys, vol. 3, 
%   pp. 17, Mar. 2015.
%   [2] https://en.wikipedia.org/wiki/Median_absolute_deviation
%   [3] https://en.wikipedia.org/wiki/Robust_measures_of_scale
%   [4] https://www.itl.nist.gov/div898/software/dataplot/refman2/auxillar/mad.htm
%
% AUTHOR:
%   Agathe BOUIS, 22/01/2024
%
% INPUTS:
%   conn_val        connected nodes' values 
%
% OUTPUTS:           
%   conn_val_filt   connected nodes' filtered out values

    
    % median 
    med = median(conn_val);
    % median absolute deviation
    MAD = median(abs(conn_val-med));
    % scaled MAD
    MAD_scaled = MAD/0.6745 ;
    % connected values z-score (median based)
    conn_z = abs((conn_val-med)./MAD_scaled);

    % Hampel identifier, values outside of +/- 3*NMAD 
    filter_val = 3;

    % node id z-score
    z_node_id = abs((node_id_val-med)./MAD_scaled);

    % min filter val
    % catch error
    % prevent filter from reaching inf when MAD=0
    if isinf(z_node_id)
        z_node_id=0;
    end
    
    % if node's z-score exceed Hampel filter, node itself is ouliter
    if z_node_id>filter_val  
        outlier_check = true;
        z_filter = filter_val;       
    else 
        outlier_check = false;
        z_filter = z_node_id;
    end
    
    % filtered values
    conn_val_filt = conn_val(conn_z<=z_filter);

end