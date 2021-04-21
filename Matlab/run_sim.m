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
sim_constants.time_step = 0.1; % s
sim_constants.tolerance = 10^-8;

% Plot formatting
plot_format.downsample_factor = 5;
plot_format.mission_name = 'LS2';
plot_format.plot_orbit_visuals = true;
plot_format.dock_plots = true;
plot_format.plot_triads = true;
% set 's', 'min', 'hours', 'days', or 'years'
plot_format.time_increments = 'hours';
plot_format = check_time_increments(plot_format);

%% Run + process sim
% sim_constants.use_euler = true; % use 312 Euler angles for attitude rather than quaternions
% sim_constants.angvel0 = deg2rad([6; 0; 0]); %rad/s, initial angular rate along principal axis
% sim_constants.simulation_time = orbital_period; % only one orbit for now
% sim_constants.perturbations_on = false; % toggle J2/drag perturbations on/off
% sim_constants.q0 = [0; sqrt(2)*0.5; 0; 0.5*sqrt(2)];
% sim_constants.q0 = A2q(OE2RTN(sim_constants.a0,sim_constants.e0,sim_constants.i0,...
%                                  sim_constants.RAAN0, sim_constants.w0, sim_constants.M0, ...
%                                  sim_constants.mu_Earth).');
sim('Propagator');

% Extract + plot data
sim_output = extract_sim_output(sim_constants, plot_format, OE, dOE_dt, ...
    w, q, e, A, ECI_positions, ECEF_positions, RTN2ECI, geod_positions);
plot_sim_output(sim_constants, sim_output, plot_format);

%% beep beep

beep;
pause(0.5);
beep;