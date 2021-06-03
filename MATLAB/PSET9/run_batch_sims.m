% close all
clear all

% Add path to subfolders
addpath('..');
addpath(fullfile('..', 'functions_and_helper_scripts'));

% Load constants
constants;

%% User input

% Simulation settings
sim_constants.simulation_time = 6*3600;
sim_constants.time_step = 0.25; % s
sim_constants.tolerance = 10^-8;

% Plot formatting
plot_format.downsample_factor = 50;
plot_format.mission_name = 'LS2';
plot_format.plot_orbit_visuals = false;
plot_format.dock_plots = true;
plot_format.plot_triads = false;
plot_format.ext_torques_verbose = false;
% set 's', 'min', 'hours', 'days', or 'years'
plot_format.time_increments = 'hours';
plot_format = check_time_increments(plot_format);

% Set batch values
plot_tag = true;
% Previously tried permutations: 
% 0.001 0.2 0.5 1 5 10
% 0.1 0.25 0.5 0.7 0.9
response_freq = [0.011]; % 0.011
damp_factor = [0.7]; % 0.7

i = length(response_freq); j=length(damp_factor);
response_freq = repmat(response_freq, [1,j]); 
damp_factor = repmat(damp_factor, [1,i])
n_sims = length(response_freq);
if n_sims ~= length(damp_factor);
    error('response_freq and damp_factor should be the same length');
end

sim_constants.state_debounce_time = 60; % sec
sim_constants.slew_large_angle_lim = deg2rad(5); % cut-off above which non-linear control is used
sim_constants.tumble_limit = deg2rad(6); % rad/s, angular velocity limit at which SC goes into detumble mode
sim_constants.target_average_window = 120; % sec, duration of moving average for target attitude

%% Run + process sims
sim_constants.attitude_perturbations_on = true;
sim_constants.orbital_perturbations_on = true;
sim_constants.sensor_noise = true;

sim_results = table;

for sim_num = 1:n_sims
    
    % Set variables
    sim_constants.k_p = sim_constants.I_princ*(response_freq(sim_num))^2;
    sim_constants.k_d = 2*damp_factor(sim_num)*sqrt(sim_constants.I_princ.*sim_constants.k_p);
    
    fprintf('\n %d. Starting sim...', sim_num);
    fprintf('\n f, d:%f, %f', response_freq(sim_num), damp_factor(sim_num));
    fprintf('\n k_p, k_d: %f\t%f', sim_constants.k_p, sim_constants.k_d);
    
    % Run it
    fprintf('\nRunning sim...');
    sim('Propagator');
    fprintf('\n\nSim complete.\n');

    % Extract + plot data if desired; extraction necessary for calculations
    sim_output = extract_sim_output(sim_constants, plot_format, ...
        OE, dOE_dt, w, w_r, dw_rdt, q, e, A, e_err, A_err,e_est, A_est, ...
        e_esterr, A_esterr, e_target, A_target, e_trueerr, A_trueerr, ...
        M_perturbations, M_drag, M_grav, M_mag, M_SRP, M_control, ...
        M_requested, m_magtor, ECI_positions, ECEF_positions, RTN2ECI, ...
        geod_positions, SC_Mode);

    if plot_tag
        fprintf('\nPlotting...\n');
        plot_format.mission_name = strcat('LS2 sim batch #', num2str(sim_num));
        plot_batch_sims(sim_constants, sim_output, plot_format);
    end

    fprintf('\nDone.\n');

end

%% beep beep, save table

for b = 1:3
    beep;
    pause(0.5);
end