close all
clear all

% Add path to subfolders
addpath(fullfile('.', 'functions_and_helper_scripts'));

% Load constants
constants;

%% User input

% Simulation settings
% Or set simulation_time in thrust_input_file.
sim_constants.simulation_time = 3.5*3600;
sim_constants.time_step = 10; % s
sim_constants.tolerance = 10^-8;

% Plot formatting
plot_format.downsample_factor = 10;
plot_format.mission_name = 'LS2';
plot_format.plot_orbit_visuals = true;
% set 's', 'min', 'hours', 'days', or 'years'
plot_format.time_increments = 'hours';
plot_format = check_time_increments(plot_format);

%% Run + process sim

sim('Propagator');

% Extract + plot data
sim_output = extract_sim_output(sim_constants, plot_format, OE, dOE_dt, ...
    ECEF_positions, geod_positions);
plot_sim_output(sim_constants, sim_output, plot_format);

%% beep beep

beep;
pause(0.5);
beep;