function [contacts_matrix, tof_matrix, pos] = constellation_contacts(simulation_param,sat_sim_type, constellation_param, isl_param, show_sat)
    % Contact Plans
    %
    % DESCRIPTION:
    %   Create contact plan (sparse matrix) for a satellite constellation. 
    %
    % INPUTS:
    %   simulation_param        {cell}
    %       start time               datetime of start of simulation
    %       stopTime                 datetime of end of simulation 
    %       sampleTime               simulation sampling/resolution   
    %   sat_sim_type            type of simulation
    %                               "Walker-Delta"
    %                               "Custom Constellation"
    %   constellation_param     {cell} of constellation parameters
    %   if "Walker-Delta"
    %       radius                  Orbital radius
    %       inclination             Inclination
    %       tot_sat                 Total Number of Satellites
    %       geometry_planes         Number of orbital planes
    %       phasing                 Phasing
    %   if "Custom"
    %       semimajoraxis           Orbital semimajor axis [length=>1]
    %       eccentriciy             Eccentricity [length=>1]
    %       inclination             Inclination [length=>1]
    %       RAAN                    RAAN [length=>1]
    %       argumentOfPeriapsis     Argument Of Periapsis [length=>1]
    %       tot_sat                 Total Number of Satellites
    %   isl_param               isl range of each satellite
    %
    % AUTHOR:
    %   Agathe BOUIS 19/04/2024

    
    % initiliase
    % number samples
    delta = simulation_param{2}-simulation_param{1};
    n_Samples = floor(seconds(delta)/simulation_param{3}) +1;
    % time lapse string
    %D = string(simulation_param{1}:seconds(simulation_param{3}):simulation_param{2}); 
    % speed of light
    c = 299792.458;
    
    %% SIMULATE
    sc = satelliteScenario(simulation_param{1},simulation_param{2},simulation_param{3});
    switch sat_sim_type
        case "Walker-Delta"
            sat = walkerDelta(sc, constellation_param{1}, constellation_param{2}, ...
                constellation_param{3}, constellation_param{4}, constellation_param{5}); ...
                n_sat = constellation_param{3};
        case "Custom Constellation"
            full_constellation_param = sat_constellation_param(constellation_param{6}, constellation_param(1:5));
            sat = satellite(sc,full_constellation_param{1},full_constellation_param{2},full_constellation_param{3}, ...
                    full_constellation_param{4},full_constellation_param{5}, full_constellation_param{6});
            n_sat = constellation_param{6};
    end
    if show_sat
        show(sat)
    end
    
    % satellite positions
    pos = states(sat,"CoordinateFrame","ecef");
    pos = pos./1000;
    
    % reframe
    pos_reframed = permute(pos, [2 3 1]);    
    pos_x =  squeeze(pos_reframed(:,:,1));
    pos_y =  squeeze(pos_reframed(:,:,2));
    pos_z =  squeeze(pos_reframed(:,:,3));
    
    % functions 
    % position difference
    diff_func = @(row) row-row.';
    % vector magnitude
    norm_func = @(row) norm(row);
    % apply to matrices
    applyToGivenRow = @(func, matrix) @(row) func(matrix(row, :));
    newApplyToRows = @(func, matrix) arrayfun(applyToGivenRow(func, matrix), 1:size(matrix,1), 'UniformOutput', false)';
    %applyToGivenMatrix = @(func, matrix_3d) @(matrix)  newApplyToRows(func, matrix_3d(:,:,matrix));
    %newApplyToMatrix = @(func, matrix_3d) arrayfun(applyToGivenMatrix(func, matrix_3d), 1:size(matrix_3d,1), 'UniformOutput', false)';
        
    % position separation between all satellites in x, y, z
    sep.X = newApplyToRows(diff_func, pos_x);
    sep.Y = newApplyToRows(diff_func, pos_y);
    sep.Z = newApplyToRows(diff_func, pos_z);
    
    % for each time step, check contact
    for i=1:n_Samples 
        % concatenate 3 dimensionsional separations back into 1 matrix
        % 3d matrix corresponds to the separation between all satellites at
        % 1 time step
        % remove lower triangle of matrix
        % matrix is symmetrical so not useful 
        % ie, x_distance sat1 -> sat2 is the same as sat2 -> sat1
        sep_mat(:,:,1) = abs(triu(cell2mat(sep.X(i))));
        sep_mat(:,:,2) = abs(triu(cell2mat(sep.Y(i))));
        sep_mat(:,:,3) = abs(triu(cell2mat(sep.Z(i))));

        % id of satellites combinations based on the reshaped matrix
        % ids not considered for matrix values of 0
        % ie, if x, y, z ===0, means that either there is no separation
        % between the two satellites (x, y, z) are the same
        % or all 0 from having lower triangle
        [id1, id2] = find(sep_mat(:,:,1)|sep_mat(:,:,2)|sep_mat(:,:,3));
    
        % reshape 3d matrix into single matrix with 3 columns x, y, z
        sep_mat_formated = reshape(sep_mat, [], 3);
        sep_mat_formated( ~any(sep_mat_formated,2), : ) = [];        

        % find magnitude of distances between sats norm(x, y, z) 
        % create sparse matrix with those magnitudes & ids
        separation = sparse(cell2mat(newApplyToRows(norm_func, sep_mat_formated)));
        separation_formated = sparse(id1,id2,separation);

        separation_formated_full = cat(1,full(separation_formated), zeros(1,n_sat));
        S = separation_formated_full+separation_formated_full.';
        isl_mat = repmat(isl_param, 1, n_sat);
        %isl_mat(logical(eye(size(isl_mat)))) =0;   
        % create contacts if separation is below isl 
        contacts_timestep_mat = S<=isl_mat;
        contacts_timestep_mat(logical(eye(size(contacts_timestep_mat)))) =0; 

        contacts_timestep = reshape(contacts_timestep_mat, [],1);
        
        % time of flight
        tof_timestep_mat = S./c;
        tof_timestep_mat = contacts_timestep_mat.*tof_timestep_mat;
        tof_timestep = reshape(tof_timestep_mat, [],1);

        % add to contact matrix
        if i == 1
            %contacts_matrix = single(contacts_current);
            contacts_matrix = contacts_timestep;
            tof_matrix = tof_timestep;
        else
            %contacts_matrix = cat(2, contacts_matrix, single(contacts_current));
            contacts_matrix = cat(2, contacts_matrix, contacts_timestep);
            tof_matrix = cat(2, tof_matrix, tof_timestep);
        end
    end
   
end