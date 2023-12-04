% Define the function to calculate fan area
function fanArea = calculateFanArea(M, gamma,T0,P0, m)
    % Constants
    R = 287.05;
    % Calculate the fan area
    fanArea = (m * sqrt(R * T0)) / (P0 * sqrt(gamma) * M * (1 + ((gamma - 1) / 2) * M^2)^(- (gamma + 1) / (2 * (gamma - 1))));

%     fanArea = (m * sqrt(R * T0)) / (P0 * sqrt(gamma) * M * ((1 + ((gamma - 1) / 2) * M^2))^(-(gamma + 1) / (2 * (gamma - 1))));
end