close all
clear all

% Add path to subfolders
addpath(fullfile('..', 'functions_and_helper_scripts'));

% Load constants
constants;

%% User input

% Simulation settings
% Or set simulation_time in thrust_input_file.
sim_constants.simulation_time = 100;
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
sim_constants.q0 = [0; 0; 0; 1];
sim_constants.angvel0 = deg2rad([0; 0; 8]);
sim_constants.use_euler = false;
sim('Propagator');

% Extract data for quaternions
sim_output_quat = extract_sim_output(sim_constants, plot_format, OE, dOE_dt, ...
    w, w_r, q, e, A, ECI_positions, ECEF_positions, RTN2ECI, geod_positions);

% Extract data for 312 Euler angles
sim_constants.use_euler = true;
sim('Propagator');
sim_output_eul = extract_sim_output(sim_constants, plot_format, OE, dOE_dt, ...
    w, w_r, q, e, A, ECI_positions, ECEF_positions, RTN2ECI, geod_positions);

%% Plot Euler angles, quaternions for each and compare error
figure('Name','Quaternions'); hold on;
time_df = downsample(sim_output_quat.time, plot_format.downsample_factor);
for i = 1:4
   subplot(2,4,i); hold on;
   title(['q_' num2str(i)]);
   xlabel('Time, hrs'); ylabel(['q_' num2str(i)]);
   plot(time_df, downsample(getfield(sim_output_quat.attitude, ['q' num2str(i)]), plot_format.downsample_factor),'-s','MarkerSize',5,'DisplayName','Quaternion');
   plot(time_df, downsample(getfield(sim_output_eul.attitude, ['q' num2str(i)]), plot_format.downsample_factor),'DisplayName','312 Euler');
   legend;
   subplot(2,4,i+4); hold on;
   title(['Error in q_' num2str(i)]);
   xlabel('Time, hrs'); ylabel(['Error in q_' num2str(i)]);
   q_err = getfield(sim_output_quat.attitude, ['q' num2str(i)])-getfield(sim_output_eul.attitude, ['q' num2str(i)]);
   plot(time_df, downsample(q_err, plot_format.downsample_factor));
end
sgtitle('Quaternion Comparison for Quaternion Param. (top) vs 312 Euler Param. (bottom)');