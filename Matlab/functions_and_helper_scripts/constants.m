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

% Satellite data
sim_constants.mass_sat = 4.93; % kg
sim_constants.effective_area_sat_max = 31.0026; % m^2
sim_constants.angvel0 = deg2rad([-6; 8; 0.1]); %rad/s, initial angular rate
sim_constants.q0 = [0; 0; 0; 1]; % principal axes initially aligned with inertial
sim_constants.I = [3.10553 -0.00011 -0.00003;
                   -0.00011 3.10289 -0.00005;
                   -0.00003 -0.00005 5.98305]; % inertia matrix, kg*m^2
% I_princ: 3x1 vector of in-order ascending principle moments of inertia, kg*m^2
% rotm: rotation matrix from principle axes to body
[sim_constants.I_princ, sim_constants.rotm] = rot_body2princ(sim_constants.I); 
sim_constants.Cdrag = 1.17; % drag coefficient
sim_constants.Cspec = 0.5; % spectral reflection coefficient for SRP
sim_constants.Cdiff = 0.5; % diffusive reflection coefficient for SRP
sim_constants.target_att_function = @LS2_target_attitude;

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
sim_constants.r_rotor = (sim_constants.rotm.')*[1;0;0]; % orientation of rotor in body coords (along body X axis)

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

% This is arbitrary.
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
sim_constants.orbital_perturbations_on = true; % add J2 and drag perturbations (true = yes)
sim_constants.attitude_perturbations_on = true; % add gravity gradient perturbations (true = yes)
