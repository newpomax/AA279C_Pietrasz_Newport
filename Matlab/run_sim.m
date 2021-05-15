close all
clear all

% Add path to subfolders
addpath(fullfile('.', 'functions_and_helper_scripts'));

% Load constants
constants;

%% User input

% Simulation settings
sim_constants.simulation_time = 0.1*3600;
sim_constants.time_step = 0.1; % s
sim_constants.tolerance = 10^-8;

% Plot formatting
plot_format.downsample_factor = 5;
plot_format.mission_name = 'LS2';
plot_format.plot_orbit_visuals = true;
plot_format.dock_plots = true;
plot_format.plot_triads = true;
plot_format.ext_torques_verbose = true;
% set 's', 'min', 'hours', 'days', or 'years'
plot_format.time_increments = 'hours';
plot_format = check_time_increments(plot_format);

%% Run + process sim
sim_constants.angvel0 = deg2rad([0; 0; 5]);
sim_constants.attitude_perturbations_on = false;
sim_constants.orbital_perturbations_on = false;
sim_constants.sensor_noise = true;
sim('Propagator');

%% Extract + plot data
sim_output = extract_sim_output(sim_constants, plot_format, OE, dOE_dt, ...
    w, w_r, q, e, A, e_err, A_err, e_target, A_target, M_perturbations, ...
    M_drag, M_grav, M_mag, M_SRP, ...
    ECI_positions, ECEF_positions, RTN2ECI, geod_positions);
plot_sim_output(sim_constants, sim_output, plot_format);

%% beep beep

beep;
pause(0.5);
beep;