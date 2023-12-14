%% CFD GEN1 PRESSURE RATIO VARIATION

% GEN1 PR vs M_flow
dataMatrix = [
    3, 1.2799, 51.5695;
    4, 1.3984, 51.5697;
    4.5, 1.4880, 51.5015;
    5, 1.5879, 51.4180;
    5.5, 1.6996, 49.9112;
    5.6, 1.7191, 48.8302;
    5.615, 1.6972, 45.0731;
    5.7, 1.2410, 0;
 
];

% Gen 1 Eff vs massflow
dataMatrix_2 = [
    3, 78.8105, 51.5695;
    4, 81.4747, 51.5697;
    4.5, 84.6385, 51.5015;
    5, 87.4143, 51.4180;
    5.5, 87.1471, 49.9112;
    5.6, 86.1781, 48.8302;
    5.615, 81.9313, 45.0731;
    5.7, 0, 0;
];

% Extracting columns 
backpressure = dataMatrix(:, 1);
fanPressureRatio = dataMatrix(:, 2);
massFlow = dataMatrix(:, 3);
eff = dataMatrix_2(:, 2);

% Plotting
figure;
plot(massFlow, fanPressureRatio, 'o-');
xlabel('Mass Flow');
ylabel('\pi');
title('Pressure Ratio');

% Plotting
figure;
plot(massFlow, eff, 'o-');
xlabel('Mass Flow');
ylabel('\eta');
title('Polytropic Efficiency');




