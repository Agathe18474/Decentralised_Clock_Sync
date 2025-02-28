function V_corrected = dV_correction(V_current, dH, rate_corrected)
    % extent of the correction of each time step 
    % Inputs: 
    % Outputs: 
    
    V_corrected = V_current - rate_corrected*dH;
    %V_corrected = V_past + rate_corrected*dH;
end