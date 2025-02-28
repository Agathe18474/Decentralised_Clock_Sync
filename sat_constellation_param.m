function sat_constellation_param_array = sat_constellation_param(totalSatellites, input_array)
    % DESCRIPTION:
    %   Create custom constellation based on Keplerian orbital parmeters.
    %   Spreads satellites equally across the number of planes
    %   corresponding to the maximum length of input parameters. 
    %   eg: if RAAN has length of 2, but all other input of size 1, will
    %   create 2 planes with different RAAN but equal inc, 
    %
    % INPUTS
    %   tot_sats                    Total number of satellites [single]
    %   input_array
    %       semimajoraxis           Orbital semimajor axis [length=>1]
    %       eccentriciy             Eccentricity [length=>1]
    %       inclination             Inclination [length=>1]
    %       RAAN                    RAAN [length=>1]
    %       argumentOfPeriapsis     Argument Of Periapsis [length=>1]
    %
    % OUTPUTS: 
    %   sat_constellation_param_array Array of parameters for all
    %                                   constellation satellites
    %
    % AUTHOR:
    %   Agathe BOUIS 19/04/2024

% shift value to avoid overlapping anomalies if satellites in the same
% orbit but different planes

% check for number of planes across which to spread the satellites 
% find the length of each param
numel_inputs = (cellfun(@numel, input_array));
numel_remainder = ones(1,numel(input_array));
numel_remainder(numel_inputs~=1) = numel_inputs(numel_inputs~=1)./max(numel_inputs);

% check for error 
% do the param of which the length is not 1 match in dimension? 
%   eg: are the non-length 1 param of equal length? 
%   if not (eg, length 3 and length 4), error
% is the number of planes (max param length) divisible by the number of
%   sats? dimensions must match to allow satellites to be distributed
%   equally across planes
if sum(numel_remainder)~=numel(input_array) || rem(totalSatellites, max(numel_inputs))~=0
    error("Ensure that number of satellites is divisible by number of planes.")
end

% create parameters for all satellites
switch_counter = totalSatellites/max(numel_inputs);
% if number of inputs not 1 - shape anomalies of satellites to avoid
% overlaps, ie so that 2 satellites in diff orbit don't have the same
% anomalies 
if max(numel_inputs)~=1
    shift = 360/((switch_counter+1)*(max(numel_inputs)+1));
    anomaly = linspace(0,360, switch_counter+1);
    anomalies = repmat(anomaly(1, 1:end-1), [max(numel_inputs) 1]) + repmat( (0:shift:shift*(max(numel_inputs)-1)).', [1 switch_counter]);
    anomalies(anomalies>360) = anomalies(anomalies>360)-360;
    trueAnomaly = reshape(anomalies.', 1, []);
else 
    anomaly = linspace(0,360, totalSatellites+1);
    trueAnomaly = anomaly(1, 1:end-1);
end

for i=1:numel(input_array)
    % if param length is 1, repeat for all satellites
    % if not, find repition pattern
    %   repeat param for satellites as per repition pattern
    if numel_inputs(i)==1
        input_array{i}=repelem(cell2mat(input_array(i)),totalSatellites);
    else
        repetion_pattern = switch_counter*ones(1,max(numel_inputs));
        input_array{i}=repelem(cell2mat(input_array(i)),repetion_pattern);
    end
end

% parameters for all constellation satellites
sat_constellation_param_array = input_array;
sat_constellation_param_array{1,6} = trueAnomaly;
end 
