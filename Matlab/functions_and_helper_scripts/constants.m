% All the constants needed for Propagator.slx

sim_constants = struct;

% Universal constants
sim_constants.M_Earth = 5.972*10^24;  % kg
sim_constants.R_Earth = 6371;         % km
sim_constants.G = 6.674 * 10^-20;     % km3 kg-1 s-2
sim_constants.mu_Earth = sim_constants.G*sim_constants.M_Earth; % km3 s-2
sim_constants.J2 = 0.00108263;        % SMAD p. 143
sim_constants.dRAAN_dt_sun_synch = 360/(365.24*24*3600); % deg s-1
sim_constants.SRP = 9.08e-6 ;         % Pa, SRP at 1AU
igrf = readmatrix('igrf13coeffs.txt');
sim_constants.igrf = igrf(2:25,28)*1e-9; % IGRF 2020 model coefficients, Tesla
% First N=4 coefficients ordered by n, then m, then g/h

% Satellite mass properties
sim_constants.mass_sat = 4.93; % kg
sim_constants.I = [3.10553 -0.00011 -0.00003;
                   -0.00011 3.10289 -0.00005;
                   -0.00003 -0.00005 5.98305]; % inertia matrix, kg*m^2
% I_princ: 3x1 vector of in-order ascending principle moments of inertia, kg*m^2
% rotm: rotation matrix from principle axes to body
[sim_constants.I_princ, sim_constants.rotm] = rot_body2princ(sim_constants.I); 

% Perturbation coefficients
sim_constants.effective_area_sat_max = 31.0026; % m^2 TODO: Replace with calculated effective area
sim_constants.Cdrag = 1.17; % drag coefficient
sim_constants.Cspec = 0.5; % spectral reflection coefficient for SRP
sim_constants.Cdiff = 0.5; % diffusive reflection coefficient for SRP
sim_constants.mag_mom_body = 4*pi*10^-7*0.1*[0.1; 0.1; 32]; % magnetic moment
sim_constants.mag_mom_princ = sim_constants.rotm*sim_constants.mag_mom_body;

% Satellite surfaces
[C, N, sim_constants.surf_areas] = get_surfacedata(); % get surface information from CSV (surface area in m^2)
sim_constants.surf_centroids = C*sim_constants.rotm; % centroid locations in principal axes, in m
sim_constants.surf_norms = N*sim_constants.rotm; % surface normals in principal axes
sim_constants.cp = sum((sim_constants.surf_areas*ones(1,size(sim_constants.surf_centroids,2)))...
                        .*sim_constants.surf_centroids)/sum(sim_constants.surf_areas); % center of pressure in princ axes

% Momentum wheel (from data sheet)
sim_constants.I_r = 0.8*0.226*(0.032^2); % kg*m^2, moment of inertia of wheel
sim_constants.w_r0 = 1; % rad/s, initial angular velocity of wheel
sim_constants.w_rmax = 0.18/sim_constants.I_r; % rad/s, maximum spin rate of wheel
sim_constants.r_rotor = (sim_constants.rotm.')*[1;0;0]; % orientation of rotor in princ coords (along body X axis)

% Magnetorquers (from data sheet)
sim_constants.m_magtor_max = 2.0; %Am^2

% Control coefficients + constants
sim_constants.k_p = 1;
sim_constants.k_d = 1;
sim_constants.state_debounce_time = 10; % sec
sim_constants.slew_lower = sind(3); % mean anomaly delta from breakpoint at which slew begins
sim_constants.slew_upper = sind(10); % mean anomaly delta from breakpoint at which slew ends
sim_constants.tumble_limit = deg2rad(15); % rad/s, angular velocity limit at which SC goes into detumble mode

% Attitude ICs
sim_constants.angvel0 = deg2rad([-6; 8; 0.1]); %rad/s, initial angular rate
sim_constants.q0 = [0; 0; 0; 1]; % principal axes initially aligned with inertial

% Attitude sensors
sim_constants.sensor_noise = true;
sim_constants.gyro_error = (deg2rad([10 10 10])/2).^2; % rad/axis
sim_constants.gyro_bias = deg2rad([0.02 0.01 -0.03]/3600); % rad/axis
sim_constants.ss_error = (deg2rad([3 3])/2).^2; % rad
sim_constants.ss_bias = deg2rad([0.7 0.2]); % rad
sim_constants.mag_error = (1E-6*[1 1 1]/2).^2; % Tesla, expect readings near 30-60E-6 Tesla
sim_constants.mag_bias = [1E-7 1E-8 -2E-8]; % Tesla
sim_constants.sensor_weights = [1 1 1]; % even weighting, three measurements made from 2 sensors (sun + magnetometer)

% EKF Input
sim_constants.Q = 1E-2*sim_constants.time_step*eye(7); % process noise
sim_constants.R = eye(8); % measurement noise 
sim_constants.R(1:3,1:3) = diag(sim_constants.gyro_error);
sim_constants.R(4:6,4:6) = diag(sim_constants.mag_error);
sim_constants.R(7:8,7:8) = diag(sim_constants.ss_error);

% Orbital elements (from
% https://secure.planetary.org/site/SPageNavigator/mission_control.html)
% Day of year	Date,       Mean Motion a           h           e
% 189.19929602  8 July 2019 14.52524232 7095.553    717.4175    0.0010951

% Epoch data 
sim_constants.start_date = [2019, 07, 08.19929602]; % [YYYY, MM, DD.dd]
sim_constants.MJD0 = mjuliandate(sim_constants.start_date); % Modified julian date at simulation start
sim_constants.vernal_equinox = [2019, 03, 20.7594907];

% Orbital elements (from data)
sim_constants.a0 = sim_constants.R_Earth + 717.4175; % km
sim_constants.e0 = 0.0010951;
sim_constants.i0 = 24;  % deg

% (arbitrarily selected)
sim_constants.RAAN0 = 0; % deg
sim_constants.w0 = 0; % deg
sim_constants.M0 = 0; % rad

sim_constants.n0 = sqrt(sim_constants.mu_Earth/sim_constants.a0^3); % rad s-1
orbital_velocity = sim_constants.n0*sim_constants.a0; % km s-1
orbital_period = 2*pi/sim_constants.n0; % s

% Default sim parameters in case they're not set elsewhere (should be set
% in run_sim or whichever script calls the sim)
sim_constants.simulation_time = .067*24*3600; % day -> hr -> s
sim_constants.time_step = 0.5; % s
sim_constants.tolerance = 10^-8;

sim_constants.use_euler = false; % use quaternions
sim_constants.use_qmethod = true; % use q-method over deterministic method
% add J2 and drag perturbations
sim_constants.orbital_perturbations_on = true;
% add gravity gradient, magnetic, drag, and SRP perturbations
sim_constants.attitude_perturbations_on = true; 
