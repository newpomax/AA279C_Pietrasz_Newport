close all
clear all

% Add path to subfolders
addpath(fullfile('.', 'functions_and_helper_scripts'));

% Load constants
constants;

%% User input

% Simulation settings
sim_constants.simulation_time = 1.5*3600;
sim_constants.time_step = 0.05; % s
sim_constants.tolerance = 10^-8;

% Plot formatting
plot_format.downsample_factor = 15;
plot_format.mission_name = 'LS2';
plot_format.plot_orbit_visuals = true;
plot_format.dock_plots = true;
plot_format.plot_triads = true;
plot_format.ext_torques_verbose = false;
% set 's', 'min', 'hours', 'days', or 'years'
plot_format.time_increments = 'hours';
plot_format = check_time_increments(plot_format);

%% Run + process sim
sim_constants.angvel0 = deg2rad([0;3;3]);
sim_constants.attitude_perturbations_on = false;
sim_constants.orbital_perturbations_on = false;
sim_constants.sensor_noise = false;
sim('Propagator');

%% Extract + plot data
sim_output = extract_sim_output(sim_constants, plot_format, OE, dOE_dt, ...
    w, w_r, q, e, A, e_err, A_err, e_est, A_est, e_esterr, A_esterr, e_target, A_target, M_perturbations, ...
    M_drag, M_grav, M_mag, M_SRP, ...
    ECI_positions, ECEF_positions, RTN2ECI, geod_positions, SC_Mode);
plot_sim_output(sim_constants, sim_output, plot_format);

%% beep beep

beep;
pause(0.5);
beep;