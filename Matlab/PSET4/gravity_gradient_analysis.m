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