function [conn_val_filt, outlier_check] = filter_Dixon(conn_val, node_val_id)
    sorted_dist = sort([conn_val; node_val_id]);
    differences = diff(sorted_dist);
    range = sorted_dist(end)-sorted_dist(1);
    Q = differences/range;
    l = length(sorted_dist);
    outlier_check = false;
    switch l
        case 3
            Q_table = 0.941;
        case 4
            Q_table = 0.765;
        case 5
            Q_table = 0.642;
        case 6
            Q_table = 0.560;
        case 7
            Q_table = 0.507;
        case 8
            Q_table = 0.468;
        case 9
            Q_table = 0.437;
        case 10
            Q_table = 0.412;
        case 11
            Q_table = 0.392;
        case 12
            Q_table = 0.376;
        case 13
            Q_table = 0.361;
        case 14
            Q_table = 0.349 ;
        case 15
            Q_table = 0.338;
        case 16
            Q_table = 0.329;
        otherwise 
            Q_table = 0.3;
    end 
    outlier_test = Q>Q_table;
    if l==2
        outlier_check = true;
        conn_val_filt = [];
    else
        if rem(l,2) ~=0
            div = (l-1)/2;
            outlier_array = zeros(1,l);
            outlier_array(1:div) = outlier_test(1:div);
            outlier_array(end-div+1:end) = outlier_test(end-div+1:end);
        else
            div = l/2;
            outlier_array = zeros(1,l);
            outlier_array(1:div+1) = outlier_test(1:div+1);
            outlier_array(end-div+1:end) = outlier_test(end-div+1:end);
        end
        
        id_ind = find(sorted_dist==node_val_id);
        if outlier_array(id_ind)==true
            outlier_check = true;
        end
        
        outlier_array(id_ind) = [];
        sorted_dist(id_ind) = [];
        conn_val_filt = sorted_dist(~outlier_array);
    end

    end
	