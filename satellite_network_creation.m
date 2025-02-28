close all 
clear all 

% Satellite Network Creation 
% Create satellite network based on user specifications 

% type of scenario 
type_scenario = "Homebrew Contacts";
%type_scenario = "Custom Constellation" ;

if type_scenario== "Homebrew Contacts"
    % if using homebrew contacts - chose type of satellite constellation
    %   custom
    %   walker delta
    %sat_sim_type = "Custom Constellation";
    sat_sim_type = "Walker-Delta";

    % simulation time parameters
    start_time = datetime(2020,6,02,8,23,0);
    stop_time = start_time + hours(1);
    sample_time = 10;

    % radius Earth
    rE = 6371;
  
    if sat_sim_type == "Walker-Delta"   
        % WALKER DELTA
        altitude = 12000;
        radius = (altitude+rE)*10^(3);
        inclination = 56;
        tot_sat = 200;
        % nb of planes
        geometryPlanes = 2;
        phasing = 1;
        sat_parameters = {radius,inclination,tot_sat,geometryPlanes,phasing};

        % Inter-Satellite Link Range
        ISL_range = 4000*ones(tot_sat, 1);
        %ISL_range = (3000*rand(n_sat,1)+ 500);

        show_sat = false;
        constellation_param = {radius, inclination, tot_sat, geometryPlanes, phasing};
            
    elseif sat_sim_type == "Custom Constellation" 
        tot_sat = 500;
        
        % CREATE RANDOM NETWORK 
        possible_inclinations = 10:1:80;
        inclination=randsample(possible_inclinations,tot_sat,true);

        possible_radius = 1000:100:2000;
        radius = randsample(possible_radius,tot_sat,true);
        semimajoraxis= (rE+radius)*10^3;

        eccentricity = 0.01*ones(1,tot_sat);
        
        possible_argumentOfPeriapsis = 0:5:360;
        argumentOfPeriapsis=randsample(possible_argumentOfPeriapsis,tot_sat,true);
        %argumentOfPeriapsis = zeros(1,tot_sat);

        possible_RAAN = 0:5:360;
        RAAN=randsample(possible_RAAN,tot_sat,true);
        %RAAN  = zeros(1,tot_sat);

        
        % CUSTOM
        %semimajoraxis= 11371e3;
        %eccentriciy = 0.01;
        %inclination = 56; 
       
        %RAAN = [0];
        % if want to go in opposite direction
        %RAAN = [0,180];
        %argumentOfPeriapsis = 0;
        show_sat = true;

        % Inter-Satellite Link Range
        ISL_range = 4000*ones(tot_sat, 1);

        constellation_param = {semimajoraxis, eccentricity, inclination, RAAN, argumentOfPeriapsis, tot_sat};
    else 
        error("Invalid Input Type")
    end
    
    simulation_param = {start_time, stop_time, sample_time};
    
    % simulate own satellite constellation 
    tic
    [A, tof, pos] = constellation_contacts(simulation_param, sat_sim_type, constellation_param, ISL_range, show_sat);
    toc
end

% test and plot newly network
A_type = class(A);
n = sqrt(size(A,1));
rdm_time_step = round(200*rand(1));
[A_current_time, tof_current_time] = time_step_adajacency(A, A_type, rdm_time_step, n, tof);

pos_plot = permute(pos, [1 3 2]); 
mode_indeg= mode(sum(A_current_time));

pos_plot_3_rings = [pos_plot(1, :, rdm_time_step); pos_plot(2, :, rdm_time_step); pos_plot(3, :, rdm_time_step)];
% filename = strcat("3_rings_positions_and_A");
% save(filename, "pos_plot_3_rings", "A_current_time)


G = digraph(A_current_time);
p = plot(G, 'XData', pos_plot(1, :, rdm_time_step), 'YData',  pos_plot(2, :, rdm_time_step), 'ZData',  pos_plot(3, :, rdm_time_step), 'MarkerSize',7, 'LineWidth', 0.1);
xlabel("X position [km]")
ylabel("Y position [km]")
zlabel("Z position [km]")
p.NodeColor  = [153, 0, 153]./256;
p.EdgeColor  = [150, 150, 150]./256;
v = [-55 -5 5];
[caz,cel] = view(v);


% view([-62.1489, 14.3582])



% % specific plots 
% view([-18,14.4])
% 
% X = ['Average In-degree on Rings is ',num2str(mode_indeg)];
% disp(X)
% 
% x = [0.66 0.57];
% y = [0.54 0.60];
% annotation('textarrow',x,y,'String', '6 Hubs')
% 
% x = [0.80 0.72];
% y = [0.35 0.40];
% annotation('textarrow',x,y,'String','12 Arcs')


% if ok, save values as 
% create string
% nb planes
if sat_sim_type == "Walker-Delta"
    nb_planes = num2str(geometryPlanes);
    indeg = num2str(mode_indeg);
    filename = strcat(nb_planes, "_rings_", indeg, "_conn_ANALYSIS");
elseif sat_sim_type== "Custom Constellation"
    nb_sats = num2str(tot_sat);
    filename = strcat(nb_sats, "_sats_random_network");
end

% add '-v7.3' extension to ensure that tof can be saved - otherwise files
% must <2GB. Sad times
%save(filename, "A", "tof", "pos_plot", '-v7.3')


