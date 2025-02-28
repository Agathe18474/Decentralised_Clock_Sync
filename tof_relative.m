function tofs_receptor = tof_relative(tof, dH, dV)
    % Time of flight (tof) in terms of H- & V-time of the receptor
    % Inputs: 
    %   tof             tof in "universal/real" time
    %   dH              receptor H-time elapsed over 1 "real" timestep
    %   dV              receptor V-time elapsed over 1 "real" timestep 
    % Outputs: 
    %   tofs_receptor   tof in terms of receptor's H- & V-times 
    
    tof_H_receptor = tof*dH;
    tof_V_receptor = tof*dV;
    
    tofs_receptor = [tof_H_receptor, tof_V_receptor];
end