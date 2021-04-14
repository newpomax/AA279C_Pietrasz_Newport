%% Compare axisymmetric numeric EulerEqs against analytical equations.
close all; clearvars; clc;

%% Constants
% Add path to subfolders
addpath(fullfile('..', 'functions_and_helper_scripts'));

% Load constants
constants;

% Higher-resolution propagation
sim_constants.time_step = 1;

%% Simulation
output_nonsymmetric = sim('EulerEqs.slx', 8000);

% Force symmetry along x & y axes
sim_constants.I_princ(2) = sim_constants.I_princ(1);

output_numeric = sim('EulerEqs.slx', 8000);
output_analytic = sim('EulerEqs_axisymmetric_analytical.slx', 8000);

% Numeric vs analytic
figure();
hold on;
plot(output_numeric.tout, output_analytic.w(:,1));
plot(output_numeric.tout, output_analytic.w(:,2));
plot(output_numeric.tout, output_analytic.w(:,3));
legend({'\omega_x', '\omega_y', '\omega_z'});
title('Euler equations analytic solution in axisymmetric case');
xlabel('time (s)');
ylabel('angular velocity (rad/s)');

% Numeric vs analytic
figure();
hold on;
plot(output_numeric.tout, output_numeric.w(:,1) - output_analytic.w(:,1));
plot(output_numeric.tout, output_numeric.w(:,2) - output_analytic.w(:,2));
plot(output_numeric.tout, output_numeric.w(:,3) - output_analytic.w(:,3));
legend({'\omega_x', '\omega_y', '\omega_z'});
title({'Comparison of Euler equations numeric integration', 'to analytic solution in axisymmetric case'});
xlabel('time (s)');
ylabel({'angular velocity error (rad/s)', '\omega_{numeric} - \omega_{analytic}'});

% Numeric axisymmetric vs nonsymmetric
figure();
hold on;
plot(output_numeric.tout, output_numeric.w(:,1) - output_nonsymmetric.w(:,1));
plot(output_numeric.tout, output_numeric.w(:,2) - output_nonsymmetric.w(:,2));
plot(output_numeric.tout, output_numeric.w(:,3) - output_nonsymmetric.w(:,3));
legend({'\omega_x', '\omega_y', '\omega_z'});
title('Comparison of Euler equations numeric integration of axisymmetric and non-symmertric cases');
xlabel('time (s)');
ylabel({'angular velocity error (rad/s)', '\omega_{axisymmetric} - \omega_{non-symmetric}'});