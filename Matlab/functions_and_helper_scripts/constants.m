% All the constants needed for Propagator.slx

sim_constants = struct;

% Universal constants
sim_constants.M_Earth = 5.972*10^24;  % kg
sim_constants.R_Earth = 6371;         % km
sim_constants.G = 6.674 * 10^-20;     % km3 kg-1 s-2
sim_constants.mu_Earth = sim_constants.G*sim_constants.M_Earth; % km3 s-2
sim_constants.J2 = 0.00108263;        % SMAD p. 143
sim_constants.dRAAN_dt_sun_synch = 360/(365.24*24*3600); % deg s-1

% Epoch data
sim_constants.start_date = [2021, 09, 16]; % [YYYY, MM, DD.dd]
sim_constants.vernal_equinox = [2021, 03, 20];

% Satellite data
sim_constants.mass_sat = 14; % kg
sim_constants.effective_area_sat_max = (1*1 + 2*6*1 + 3*1)*.01; % m^2
sim_constants.effective_area_sat_min = (2*3)*.01;
sim_constants.thruster_total_impulse = 12000; % N s

% Orbital elements
sim_constants.a0 = sim_constants.R_Earth + 550; % km
sim_constants.e0 = 0;
sim_constants.i0 = 97.58;  % deg
sim_constants.MLT0 = 360-22.5; % deg -- 10:30am @ vernal equinox
sim_constants.RAAN_offset0 = sim_constants.dRAAN_dt_sun_synch*seconds( ...
        datetime(sim_constants.start_date) - ...
            datetime(sim_constants.vernal_equinox) ...
        ); % deg -- difference between RAAN0 and MLT0.
sim_constants.RAAN0 = rem(sim_constants.MLT0 + ...
    sim_constants.RAAN_offset0, 360);   

sim_constants.w0 = 359.99; % deg
% As long as we're assuming e ~ 0, w0 doesn't matter. 
% And so we'll set it to 0, which means that if dw/dt = 0,
% the argument of latitude u = nu. But, because dw/dt from J2 < 0, we'll
% actually set it to ~360 so that the plot doesn't need to wrap from 0 to
% 360 at the beginning and we can better visualize the changes.

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
