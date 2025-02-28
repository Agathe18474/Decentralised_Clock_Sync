function [A, tof, n, sat_list_names, names, conn_start,test_step_size, conn_max, pos_plot] = ring_setup(nb_rings)
    % Ring Setup - initialise ring topology
    % Inputs:   
    %   nb_rings                # rings 
    % Outputs: 
    %   A                       adjacency matrices of network for simulation 
    %   tof                     tof matrices of network for simulation
    %   sat_list_names          name of networks, nb_ring_conn   
    %   names                   connectivity on rings
    %   conn_test               array of conn for each ring

    switch nb_rings
        case 1 
            conn_start= 2; 
            test_step_size = 2;
            conn_max = 8;
            %names = {'Conn 2', 'Conn 4', 'Conn 10','Conn 20', 'Conn 30', 'Conn 40'};
            %names = {'Conn 2', 'Conn 4','Conn 6','Conn 8','Conn 10','Conn 12','Conn 14','Conn 16','Conn 18'};
            names = {'Conn 2', 'Conn 4','Conn 6','Conn 8'};
            n = 100;

            % initialise for first ring
            altitude = 2000;
            rE = 6371;
            r = (altitude+rE)*10^(3);
            k = conn_start;
            [A, tof] = ring_graph(r, n, k);
            sat_list_names = [];  

            % half ring
            %[A, tof] = half_ring_graph(r, n, k);
            % G = digraph(A)
            % plot(G)
            % arc
             % [A, tof] = half_ring_graph(r, n*2, k);
            %[A, tof] = ring_graph(r, n*2, k);
            % A = A(1:round(n), 1:round(n));
            % tof = tof(1:round(n), 1:round(n));
            %avg_tof = sum(tof, 1)./sum(tof~=0,1); 

        case 2
            %data = "2_rings_2_conn_ANALYSIS.mat";
            data = "2_rings_4_conn_ANALYSIS.mat";
            sat_data=load(data);
            A = sat_data.A;
            tof = sat_data.tof;
            pos_plot = sat_data.pos_plot;
            n = sqrt(size(A,1));
            %tof = 0016*ones(n);
            conn_start= 2; 
            test_step_size=2;
            conn_max = 8;
            % sat_list_names = ["2_rings_2_conn_ANALYSIS","2_rings_4_conn_ANALYSIS","2_rings_6_conn_ANALYSIS","2_rings_8_conn_ANALYSIS"];%,"2_rings_10_conn_ANALYSIS" ];
            % names = {'Conn 2', 'Conn 4','Conn 6','Conn 8','Conn 10'};
            sat_list_names = ["2_rings_4_conn_ANALYSIS","2_rings_6_conn_ANALYSIS","2_rings_8_conn_ANALYSIS"];%,"2_rings_10_conn_ANALYSIS" ];
            names = {'Conn 4','Conn 6','Conn 8','Conn 10'};
            
            
            % data = "2_rings_4_conn.mat";
            % sat_data=load(data);            
            % A = sat_data.A_temp;
            % n = sqrt(size(A,1));
            % tof = zeros(n,length(A));
            % pos_plot = 0;
            % conn_start= 4; 
            % test_step_size=2;
            % conn_max = 8;
            % sat_list_names = ["2_rings_4_conn","2_rings_6_conn","2_rings_8_conn"];
            % 
            % names = {'Conn 4','Conn 6','Conn 8','Conn 10'};
        case 3
            data = "3_rings_2_conn_ANALYSIS.mat";
            sat_data=load(data);
            A = sat_data.A;
            tof = sat_data.tof;
            pos_plot = sat_data.pos_plot;
            n = sqrt(size(A,1));
            %tof = 0.0016*ones(n);
            conn_start= 2; 
            test_step_size=2;
            conn_max = 6;
            sat_list_names = ["3_rings_2_conn_ANALYSIS","3_rings_4_conn_ANALYSIS","3_rings_6_conn_ANALYSIS"];%,"3_rings_8_conn_ANALYSIS" ];
            names = {'Conn 2', 'Conn 4','Conn 6'};%,'Conn 8'};
        otherwise  
            error("Number of Rings must be 1, 2, or 3.")
    end
end