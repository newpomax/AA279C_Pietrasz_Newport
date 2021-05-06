%% Perturbation modeling
close all
clear all

% Add path to subfolders
addpath(fullfile('..', 'functions_and_helper_scripts'));

% Load constants
constants;

%% User input

% Simulation settings
sim_constants.simulation_time = 0.5*3600;
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
sim_constants.q0 = [0; sqrt(2)/2; 0; sqrt(2)/2]; % point z in direction of motion, more or less
sim_constants.angvel0 = deg2rad([0; 0; 8]); % spin about z
sim('Propagator');

% Extract + plot data
sim_output = extract_sim_output(sim_constants, plot_format, OE, dOE_dt, ...
    w, w_r, q, e, A,e_err, A_err, e_target, A_target, M_perturbations, ...
    M_drag, M_grav, M_mag, M_SRP,ECI_positions, ECEF_positions, RTN2ECI, geod_positions);
plot_sim_output(sim_constants, sim_output, plot_format);

time_label = ['Time [', num2str(plot_format.time_increments), ']'];
df = plot_format.downsample_factor;
time_df = downsample(sim_output.time,df);
%% Plotting
figure('Name','Torque Comparison'); hold on;
labs = {'x', 'y', 'z'};
for i = 1:3
    subplot(1,4,i); 
    semilogy(time_df,abs(downsample(M_grav.Data(:,i),df)),'DisplayName','Gravitational');
    hold on;
    semilogy(time_df,abs(downsample(M_drag.Data(:,i),df)),'DisplayName','Drag');
     semilogy(time_df,abs(downsample(M_mag.Data(:,i),df)),'DisplayName','Magnetic Torque');
    semilogy(time_df,abs(downsample(M_SRP.Data(:,i),df)),'DisplayName','SRP');
    xlabel(time_label); ylabel('Torque, N\cdot m');
    title_text = ['M_' labs{i} ' of ', plot_format.mission_name, ' in s/c principle axes'];
    title(title_text); legend;
end
subplot(1,4,4); 
semilogy(time_df,downsample(sqrt(sum(M_grav.Data.^2, 2)), df),'DisplayName','Gravitational');
hold on;
semilogy(time_df,downsample(sqrt(sum(M_drag.Data.^2, 2)), df),'DisplayName','Drag');
semilogy(time_df,downsample(sqrt(sum(M_mag.Data.^2, 2)), df),'DisplayName','Magnetic Torque');
semilogy(time_df,downsample(sqrt(sum(M_SRP.Data.^2, 2)), df),'DisplayName','SRP');
xlabel(time_label); ylabel('Torque, N\cdot m');
title_text = ['|M| of ', plot_format.mission_name, ' in s/c principle axes'];
title(title_text); legend('location','best');

% Predicted values from Wertz
cp = norm(sim_constants.cp); % center of pressure distance from CoM
rho = atmospheric_model(sim_constants.R_Earth, ECI_positions.Data(1,:));
theta_max = pi/3;
Mg_max = 1.5*sim_constants.mu_Earth/(norm(ECI_positions.Data(1,:))^3)...
            *(sim_constants.I_princ(3)-sim_constants.I_princ(1))*sin(2*theta_max);
Msrp_max = cp*sim_constants.SRP*sim_constants.effective_area_sat_max*1.6 ;
Mdrag_max = 0.5*cp*rho*sim_constants.Cdrag*sim_constants.effective_area_sat_max*(norm(1000*ECI_velocities.Data(1,:))^2);
fprintf('Wertz Grav Torque: %0.3e Nm \n Wertz SRP Torque: %0.3e Nm \n Wertz Drag Torque: %0.3e Nm \n', Mg_max, Msrp_max, Mdrag_max);
%% beep beep

beep;
pause(0.5);
beep;

% atmo model
function rho = atmospheric_model(R_Earth, rECI)
    a = norm(rECI); % distance from Earth center in km

    % Atmospheric density
    % From http://www.braeunig.us/space/atmos.htm
    % Density @ reference altitude h0 = 550km
    rho_0 = (2.14*10^-13); % kg m-3 -> kg km-3, 
    h0 = 550 + R_Earth; % km, reference altitude -> radius
    H = 68.7; % km, atmospheric scale height
    
    rho = rho_0*exp(-(a-h0)/H); % kg km-3
    rho = rho; % kg m-3
end