%% Perturbation modeling
close all
clear all

% Add path to subfolders
addpath(fullfile('..', 'functions_and_helper_scripts'));

% Load constants
constants;

%% User input

% Simulation settings
sim_constants.simulation_time = 0.1*3600;
sim_constants.time_step = 0.01; % s
sim_constants.tolerance = 10^-8;

% Plot formatting
plot_format.downsample_factor = 10;
plot_format.mission_name = 'LS2';
plot_format.plot_orbit_visuals = true;
plot_format.dock_plots = true;
plot_format.plot_triads = true;
% set 's', 'min', 'hours', 'days', or 'years'
plot_format.time_increments = 'hours';
plot_format = check_time_increments(plot_format);

%% Run + process sim
sim_constants.angvel0 = deg2rad([0; 0; 1]);
sim_constants.attitude_perturbations_on = false;
sim_constants.orbital_perturbations_on = false;
sim('Propagator');

% Extract + plot data
sim_output = extract_sim_output(sim_constants, plot_format, OE, dOE_dt, ...
    w, w_r, q, e, A, e_err, A_err, e_target, A_target, M_perturbations, ...
    ECI_positions, ECEF_positions, RTN2ECI, geod_positions);
plot_sim_output(sim_constants, sim_output, plot_format);

time_label = ['Time [', num2str(plot_format.time_increments), ']'];
df = plot_format.downsample_factor;
time_df = downsample(sim_output.time,df);

% Plotting
figure('Name','Attitude Error'); hold on;
time_label = ['Time [', num2str(plot_format.time_increments), ']'];
df = plot_format.downsample_factor;
time_df = downsample(sim_output.time,df);
labs = {'\phi', '\theta', '\psi'};
for i = 1:3
   subplot(1,3,i); hold on;
   plot(time_df,wrapTo180(rad2deg(downsample(e_target.Data(:,i),df))),'MarkerSize',5,'DisplayName','Target Angle');
   plot(time_df,wrapTo180(rad2deg(downsample(e.Data(:,i),df))),'DisplayName','Attitude Angle');
   plot(time_df,wrapTo180(rad2deg(downsample(e_err.Data(:,i),df))),'DisplayName','Error Angle');
   legend; ylabel([labs{i} ', deg']); xlabel(time_label);
   title_text = [labs{i} ' of ', plot_format.mission_name, ' for various DCMs'];
   title(title_text); legend('location','best');
   if i>1
      ylim([-5 5]); 
   end
end

%% beep beep

beep;
pause(0.5);
beep;