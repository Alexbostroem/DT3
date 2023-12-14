frac_ht = [0.000; 0.121; 0.233; 0.339; 0.439; 0.537; 0.632; 0.726; 0.818; 0.909; 1.000];


%% Inlet


%% IGV


%% Rotor
blade_angle_rotor_te = [43.69; 47.46; 50.82; 52.72; 54.41; 55.83; 56.95; 57.78; 58.31; 58.13; 57.20];

% Create plots
figure;

% Plot 1: 
subplot(2, 2, 1);
plot(blade_angle_rotor_te, frac_ht, '-o');
title('Outflow angle of rotor (relative)');
xlabel('deg');
ylabel('Fraction of blade height');

%% Stator

%% Outlet