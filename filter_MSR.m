function [conn_val_filt] = filter_MSR(conn_val, i_dis, node_id_val)
% Synchronous Mean Subsequence Reduced (MSR) Filter
%
% DESCRIPTION: 
%   MSR algorithm: resilient consensus strategy against an arbitrary 
%   nb of faults (F) given a specified robustness (F+1,F+1).
% REFERENCES:
% [1] LeBlanc, H et al. "Resilient asymptotic consensus in robust networks." 
% IEEE Journal on Selected Areas in Communications (2013): 766-781.
%
% AUTHOR:
%   Agathe BOUIS, 22/01/2024
%
% INPUTS:
%   conn_val        connected nodes' values 
%   i_dis           array of disruptor node ids
%   node_id_val     val of node doing the filtering
%
% OUTPUTS:           
%   conn_val_filt   connected nodes' filtered out values

F = length(i_dis);
conn_val = sort(conn_val);

nb_val_smaller_id = length(conn_val(conn_val<node_id_val));

% remove F smallest values or remove "counter" if there are less
% than F values smaller than x_i
F_min = min([nb_val_smaller_id F]);

nb_val_bigger_id = length(conn_val(conn_val>node_id_val));

% remove F largest values or remove "counter" if there are less
% than F values larger than x_i
F_max = min([nb_val_bigger_id F]);

% remove from array F_min and F_max val array
conn_val_filt = conn_val(1,(F_min+1):(end-F_max));

end