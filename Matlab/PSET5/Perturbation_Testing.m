%% Perturbation modeling
close all
clear all

% Add path to subfolders
addpath(fullfile('..', 'functions_and_helper_scripts'));

% Load constants
constants;

%% User input

% Simulation settings
sim_constants.simulation_time = 1.5*3600;
sim_constants.time_step = 0.1; % s
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
sim_constants.angvel0 = deg2rad([0; 0; 8]);
sim('Propagator');

% Extract + plot data
sim_output = extract_sim_output(sim_constants, plot_format, OE, dOE_dt, ...
    w, w_r, q, e, A,e_err, A_err, e_target, A_target, M_perturbations, ...
    ECI_positions, ECEF_positions, RTN2ECI, geod_positions);
plot_sim_output(sim_constants, sim_output, plot_format);

time_label = ['Time [', num2str(plot_format.time_increments), ']'];
df = plot_format.downsample_factor;
time_df = downsample(sim_output.time,df);
% Plotting
figure('Name','Torque Comparison'); hold on;
labs = {'x', 'y', 'z'};
for i = 1:3
    subplot(1,4,i); 
    semilogy(time_df,abs(downsample(M_grav.Data(:,i),df)),'DisplayName','Gravitational');
    hold on;
    semilogy(time_df,abs(downsample(M_drag.Data(:,i),df)),'DisplayName','Drag');
%     semilogy(time_df,abs(downsample(M_mag.Data(:,i),df)),'DisplayName','Magnetic Torque');
    semilogy(time_df,abs(downsample(M_SRP.Data(:,i),df)),'DisplayName','SRP');
    xlabel(time_label); ylabel('Torque, N\cdot m');
    title_text = ['M_' labs{i} ' of ', plot_format.mission_name, ' in s/c principle axes'];
    title(title_text); legend;
end
subplot(1,4,4); 
semilogy(time_df,downsample(sqrt(sum(M_grav.Data.^2, 2)), df),'DisplayName','Gravitational');
hold on;
semilogy(time_df,downsample(sqrt(sum(M_drag.Data.^2, 2)), df),'DisplayName','Drag');
%     semilogy(time_df,downsample(sqrt(sum(M_mag.Data.^2, 2)), df),'DisplayName','Magnetic Torque');
semilogy(time_df,downsample(sqrt(sum(M_SRP.Data.^2, 2)), df),'DisplayName','SRP');
xlabel(time_label); ylabel('Torque, N\cdot m');
title_text = ['|M| of ', plot_format.mission_name, ' in s/c principle axes'];
title(title_text); legend('location','best');


%% beep beep

beep;
pause(0.5);
beep;