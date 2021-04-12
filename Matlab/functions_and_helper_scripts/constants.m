% All the constants needed for Propagator.slx

sim_constants = struct;

% Universal constants
sim_constants.M_Earth = 5.972*10^24;  % kg
sim_constants.R_Earth = 6371;         % km
sim_constants.G = 6.674 * 10^-20;     % km3 kg-1 s-2
sim_constants.mu_Earth = sim_constants.G*sim_constants.M_Earth; % km3 s-2
sim_constants.J2 = 0.00108263;        % SMAD p. 143
sim_constants.dRAAN_dt_sun_synch = 360/(365.24*24*3600); % deg s-1

% Satellite data
sim_constants.mass_sat = 4.93; % kg
sim_constants.effective_area_sat_max = 31.0026; % m^2
sim_constants.angvel0 = deg2rad([-6; 8; 0.1]); %rad/s, initial angular rate
sim_constants.I = [3.10553 -0.00011 -0.00003;
                   -0.00011 3.10289 -0.00005;
                   -0.00003 -0.00005 5.98305]; % inertia matrix, kg*m^2
[sim_constants.I_princ, sim_constants.rotm] = rot_body2princ(sim_constants.I); 
% I_princ: 3x1 vector of in-order ascending principle moments of inertia, kg*m^2
% rotm: rotation matrix from principle axes to body

% Orbital elements (from
% https://secure.planetary.org/site/SPageNavigator/mission_control.html)
% Day of year	Date,       Mean Motion a           h           e
% 189.19929602  8 July 2019 14.52524232 7095.553    717.4175    0.0010951

% Epoch data (Currently unused)
sim_constants.start_date = [2019, 07, 08.19929602]; % [YYYY, MM, DD.dd]
sim_constants.vernal_equinox = [2019, 03, 20.7594907];

% Orbital elements (from data)
sim_constants.a0 = sim_constants.R_Earth + 717.4175; % km
sim_constants.e0 = 0.0010951;
sim_constants.i0 = 24;  % deg

% (arbitrarily selected) % TODO
sim_constants.RAAN0 = 0; % deg %
sim_constants.w0 = 0; % deg

% This is arbitrary.
sim_constants.M0 = 0; % rad
sim_constants.n0 = sqrt(sim_constants.mu_Earth/sim_constants.a0^3); % rad s-1
orbital_velocity = sim_constants.n0*sim_constants.a0; % km s-1
orbital_period = 2*pi/sim_constants.n0; % s

% Default sim parameters in case they're not set elsewhere (should be set
% in run_sim or the thrust_input_file)
sim_constants.simulation_time = .067*24*3600; % day -> hr -> s
sim_constants.time_step = 10; % s
sim_constants.tolerance = 10^-8;
