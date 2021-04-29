close all
clear all

% Add path to subfolders
addpath(fullfile('..', 'functions_and_helper_scripts'));

% Load constants
constants;

%% User input

% Simulation settings
% Or set simulation_time in thrust_input_file.
sim_constants.simulation_time = 2*3600;
sim_constants.time_step = 0.05; % s
sim_constants.tolerance = 10^-8;

% Plot formatting
plot_format.downsample_factor = 2;
plot_format.mission_name = 'LS2';
plot_format.plot_orbit_visuals = true;
plot_format.dock_plots = true;
plot_format.plot_triads = true;
% set 's', 'min', 'hours', 'days', or 'years'
plot_format.time_increments = 'hours';
plot_format = check_time_increments(plot_format);

%% Select test
run_test = 'Other_Stability';

sim_constants.w_r0 = 1; % rad/s, very slow spin for testing
switch run_test
    case 'Eq_X'
        sim_constants.r_rotor = [1; 0; 0];
        sim_constants.angvel0 = deg2rad([7; 0; 0]);
    case 'Eq_Y'
        sim_constants.r_rotor = [0; 1; 0];
        sim_constants.angvel0 = deg2rad([0; 7; 0]);
    case 'Eq_Z'
        sim_constants.r_rotor = [0; 0; 1];
        sim_constants.angvel0 = deg2rad([0; 0; 7]);
    case 'RTN'
        sim_constants.perturbations_on = false; % Turn off J2
        sim_constants.r_rotor = [0; 0; 1];
        sim_constants.angvel0 = deg2rad([0; 0; 7]); % Along N
        sim_constants.q0 = A2q(OE2RTN(sim_constants.a0,sim_constants.e0,sim_constants.i0,...
                sim_constants.RAAN0,sim_constants.w0,sim_constants.M0,sim_constants.mu_Earth).'); 
                % initially aligned with RTN
    case 'Stab_X'
        sim_constants.r_rotor = [1; 0; 0];
        sim_constants.angvel0 = deg2rad([7; 0.01; 0.01]);
    case 'Stab_Y'
        sim_constants.r_rotor = [0; 1; 0];
        sim_constants.angvel0 = deg2rad([0.01; 7; 0.01]);
    case 'Stab_Z'
        sim_constants.r_rotor = [0; 0; 1];
        sim_constants.angvel0 = deg2rad([0.1; 0.1; 7]);
    case 'Int_Stability'
        sim_constants.r_rotor = [0; 1; 0];
        sim_constants.angvel0 = deg2rad([0.01;7;0.01]); %rad/s, along intermediate axis
        I = sim_constants.I_princ; w0 = sim_constants.angvel0;
        sim_constants.w_r0 = 10*(I(3)-I(2))*w0(2)/sim_constants.I_r; % ensure stability
        sim_constants.w_rmax = abs(sim_constants.w_r0);
    case 'Other_Stability'
        sim_constants.r_rotor = [0; 0; 1];
        sim_constants.angvel0 = deg2rad([0.01;0.01;0.01]); %rad/s, along z axis w/o spin
        I = sim_constants.I_princ; w0 = sim_constants.angvel0;
        sim_constants.w_r0 = 100; % ensure stability ( anything other than 0 is sufficient )
        sim_constants.w_rmax = abs(sim_constants.w_r0);
end



%% Run + process sim
sim('Propagator');

% Extract + plot data
sim_output = extract_sim_output(sim_constants, plot_format, OE, dOE_dt, ...
    w, w_r, q, e, A, ECI_positions, ECEF_positions, RTN2ECI, geod_positions);
plot_sim_output(sim_constants, sim_output, plot_format);

%% beep beep

beep;
pause(0.5);
beep;