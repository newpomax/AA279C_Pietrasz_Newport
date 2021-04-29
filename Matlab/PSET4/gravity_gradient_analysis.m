close all
clear all

% Add path to subfolders
addpath('..');
addpath(fullfile('..', 'functions_and_helper_scripts'));

% Load constants
constants;

%% User input

% Simulation settings
sim_constants.simulation_time = 2*3600; 
sim_constants.time_step = 0.05; % s
sim_constants.tolerance = 10^-8;

% Turn rotor off
sim_constants.w_r0 = 0;

% Plot formatting
plot_format.mission_name = 'LS2';
plot_format.downsample_factor = 2;
plot_format.plot_orbit_visuals = false;
plot_format.dock_plots = true;
plot_format.plot_triads = true;
% set 's', 'min', 'hours', 'days', or 'years'
plot_format.time_increments = 'hours';
plot_format = check_time_increments(plot_format);

%% Toggles for different cases
sim_constants.orbital_perturbations_on = false;
sim_constants.attitude_perturbations_on = true;

% Stable config -- align principal axes with RTN: 
% ECI_A_RTN = OE2RTN(sim_constants.a0, sim_constants.e0, sim_constants.i0, ...
%     sim_constants.RAAN0, sim_constants.w0, sim_constants.M0, ...
%     sim_constants.mu_Earth)';
% sim_constants.q0 = A2q(ECI_A_RTN);
% sim_constants.angvel0 = [0; 0; sim_constants.n0]; %rad/s

%% Quick verification calcs

format longg

I = sim_constants.I_princ;
R = [7080.6, 0, 0]'; % In this calc, principal aligned with RTN

M_est_aligned = (3*sim_constants.mu_Earth/norm(R)^5)* ...
    [(I(3) - I(2))*R(2)*R(3);
     (I(1) - I(3))*R(3)*R(1);
     (I(2) - I(1))*R(1)*R(2)]

R = rotz(45)*roty(45)*rotx(45)*R; % Rotate to arbitrary config
M_est_nonaligned = (3*sim_constants.mu_Earth/norm(R)^5)* ...
    [(I(3) - I(2))*R(2)*R(3);
     (I(1) - I(3))*R(3)*R(1);
     (I(2) - I(1))*R(1)*R(2)]
 
 %% Stability
 
kR = (I(3)-I(2))/I(1)
kT = (I(3)-I(1))/I(2)
kN = (I(2)-I(1))/I(3)

figure();
hold on;

% unstable pitch
patch([-1, -1, 1], [-1, 1, 1], 'y', 'FaceAlpha', .4);
text(0.1, 0.7, 'Unstable pitch'); 

% unstable yaw
patch([0, 0, 1, 1], [-1, 0, 0, -1], 'b', 'FaceAlpha', .4);
patch([0, 0, -1, -1], [1, 0, 0, 1], 'b', 'FaceAlpha', .4);
text(0.4, -0.5, 'Unstable yaw'); 

% unstable roll
fp = fimplicit(@(x,y) 1 + 3*x + x.*y - 4*sqrt(x.*y), [-1, 0, -1, 0], 'k'); 
fp_x = fp.XData;
fp_y = fp.YData;
patch([fp_x fliplr(fp_x)], [fp_y -1*ones(1,length(fp_y))], ...
    'r', 'FaceAlpha', .4, 'EdgeAlpha', 0)
% Estimate the rest of the curve with a parabola
unstable_roll_x = [-1:0.01:-(1/3) -(1/3)];
unstable_roll_y = -0.2*(unstable_roll_x+1/3).^2;
plot(unstable_roll_x, unstable_roll_y, 'k');
patch([unstable_roll_x, fliplr(unstable_roll_x)], ...
    [unstable_roll_y -1*ones(1,length(unstable_roll_y))], ...
    'r', 'FaceAlpha', .4, 'EdgeAlpha', 0);
text(-0.5, -0.75, 'Unstable roll'); 

% Satellite
plot([kT], [kR], 'k.', 'MarkerSize', 25);
text_string = strcat('LS2 (-', num2str(round(kT, 5)), ', ', ...
    num2str(round(kR, 5)), ')');
text(0.6, 0.93, text_string); 

% Alternative options
% Still need to finish

title('Gravity gradient stability diagram');
xlabel('k_T');
ylabel('k_R');
xlim([-1,1]);
ylim([-1,1]);

%% Run + process sim
sim('Propagator');

% Extract + plot data
sim_output = extract_sim_output(sim_constants, plot_format, OE, dOE_dt, ...
    w, w_r, q, e, A, M_perturbations, ECI_positions, ECEF_positions, ...
    RTN2ECI, geod_positions);
plot_sim_output(sim_constants, sim_output, plot_format);

%% beep beep

beep;
pause(0.5);
beep;