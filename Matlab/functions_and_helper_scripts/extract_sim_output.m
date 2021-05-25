function sim_output = extract_sim_output(sim_constants, plot_format, ...
    OE, dOE_dt, w, w_r, dw_rdt, q, e, A, e_err, A_err,e_est, A_est, e_esterr, A_esterr, e_targ, A_targ, ...
    e_trueerr,A_trueerr,M_perturbations, M_drag, M_grav, M_mag, M_SRP, M_control,M_requested,m_magtor, ...
    ECI_positions, ECEF_positions, RTN2ECI, geod_positions, SC_Mode)
    %% Post-process data output from simulation for plotting and calcs

    sim_output = struct;
    sim_output.time = OE.time/plot_format.seconds_to_increment;

    %% Postions for orbit visualizations
    sim_output.positions.x_123 = ECI_positions.Data(:,1);
    sim_output.positions.y_123 = ECI_positions.Data(:,2);
    sim_output.positions.z_123 = ECI_positions.Data(:,3);
    
    sim_output.positions.x_XYZ = ECEF_positions.Data(:,1);
    sim_output.positions.y_XYZ = ECEF_positions.Data(:,2);
    sim_output.positions.z_XYZ = ECEF_positions.Data(:,3);

    sim_output.positions.latitude_d = geod_positions.Data(:,1);
    sim_output.positions.longitude_d = geod_positions.Data(:,2);
    
    sim_output.positions.RTN2ECI = permute(RTN2ECI.Data, [3 1 2]);

    %% For OE & dOE/dt plots

    sim_output.OE.a = OE.Data(:,1); % km
    sim_output.OE.e = OE.Data(:,2); % deg
    sim_output.OE.i = OE.Data(:,3); % deg

    sim_output.OE.RAAN = OE.Data(:,4); % deg
    sim_output.OE.w = OE.Data(:,5); % deg
    sim_output.OE.E = OE.Data(:,6); % rad

    sim_output.dOE.a = dOE_dt.Data(:,1); % km
    sim_output.dOE.e = dOE_dt.Data(:,2); % deg s-1
    sim_output.dOE.i = dOE_dt.Data(:,3); % deg s-1

    sim_output.dOE.RAAN = dOE_dt.Data(:,4); % deg s-1
    sim_output.dOE.w = dOE_dt.Data(:,5); % deg s-1
    sim_output.dOE.E = dOE_dt.Data(:,6); % deg s-1
    
    %% For momentum wheel & other control
    
    sim_output.attitude.w_r = w_r.Data; % rad/s, angular velocity of momentum wheel
    sim_output.attitude.dw_rdt = dw_rdt.Data; %rad/s^2, ang acceleration of momentum wheel
    sim_output.attitude.L_r = sim_constants.I_r*sim_output.attitude.w_r; % Nm*s, angular momentum of momentum wheel
    
    sim_output.attitude.m_magtor = m_magtor.Data; % Am, dipole moment applied by mag torquer
    
    sim_output.attitude.M_req = M_requested.Data; % Torque requested by ADCS, not actuator-limited
    sim_output.attitude.M_c = M_control.Data; % True applied control torque
    
    %% For environmental perturbations
    
    sim_output.ext_torques.Mx = M_perturbations.Data(:,1); % rad s-1, principal axes
    sim_output.ext_torques.My = M_perturbations.Data(:,2); % rad s-1, principal axes
    sim_output.ext_torques.Mz = M_perturbations.Data(:,3); % rad s-1, principal axes
    
    sim_output.ext_torques.drag_Mx = M_drag.Data(:,1); % rad s-1, principal axes
    sim_output.ext_torques.drag_My = M_drag.Data(:,2); % rad s-1, principal axes
    sim_output.ext_torques.drag_Mz = M_drag.Data(:,3); % rad s-1, principal axes
    
    sim_output.ext_torques.grav_Mx = M_grav.Data(:,1); % rad s-1, principal axes
    sim_output.ext_torques.grav_My = M_grav.Data(:,2); % rad s-1, principal axes
    sim_output.ext_torques.grav_Mz = M_grav.Data(:,3); % rad s-1, principal axes
    
    sim_output.ext_torques.SRP_Mx = M_SRP.Data(:,1); % rad s-1, principal axes
    sim_output.ext_torques.SRP_My = M_SRP.Data(:,2); % rad s-1, principal axes
    sim_output.ext_torques.SRP_Mz = M_SRP.Data(:,3); % rad s-1, principal axes
    
    sim_output.ext_torques.mag_Mx = M_mag.Data(:,1); % rad s-1, principal axes
    sim_output.ext_torques.mag_My = M_mag.Data(:,2); % rad s-1, principal axes
    sim_output.ext_torques.mag_Mz = M_mag.Data(:,3); % rad s-1, principal axes
    
    %% For attitude plots in principal axes
    sim_output.attitude.SC_Mode = SC_Mode.Data;
    sim_output.attitude.wx = w.Data(:,1); % rad s-1, principal axes
    sim_output.attitude.wy = w.Data(:,2); % rad s-1, principal axes
    sim_output.attitude.wz = w.Data(:,3); % rad s-1, principal axes
    
    sim_output.attitude.Lx = sim_constants.I_princ(1)*sim_output.attitude.wx + sim_output.attitude.L_r*sim_constants.r_rotor(1); % kg m^2/s, principal axes
    sim_output.attitude.Ly = sim_constants.I_princ(2)*sim_output.attitude.wy + sim_output.attitude.L_r*sim_constants.r_rotor(2); % kg m^2/s, principal axes
    sim_output.attitude.Lz = sim_constants.I_princ(3)*sim_output.attitude.wz + sim_output.attitude.L_r*sim_constants.r_rotor(3); % kg m^2/s, principal axes
    
    sim_output.attitude.q1 = q.Data(:,1); % unitless
    sim_output.attitude.q2 = q.Data(:,2); % unitless
    sim_output.attitude.q3 = q.Data(:,3); % unitless
    sim_output.attitude.q4 = q.Data(:,4); % unitless
    
    sim_output.attitude.phi = rad2deg(e.Data(:,1)); % deg
    sim_output.attitude.theta = rad2deg(e.Data(:,2)); % deg
    sim_output.attitude.psi = rad2deg(e.Data(:,3)); % deg
    
    sim_output.attitude.err_phi = rad2deg(e_err.Data(:,1)); % deg
    sim_output.attitude.err_theta = rad2deg(e_err.Data(:,2)); % deg
    sim_output.attitude.err_psi = rad2deg(e_err.Data(:,3)); % deg
    
    sim_output.attitude.est_phi = rad2deg(e_est.Data(:,1)); % deg
    sim_output.attitude.est_theta = rad2deg(e_est.Data(:,2)); % deg
    sim_output.attitude.est_psi = rad2deg(e_est.Data(:,3)); % deg
    
    sim_output.attitude.esterr_phi = rad2deg(e_esterr.Data(:,1)); % deg
    sim_output.attitude.esterr_theta = rad2deg(e_esterr.Data(:,2)); % deg
    sim_output.attitude.esterr_psi = rad2deg(e_esterr.Data(:,3)); % deg
    
    e_targ.Data = squeeze(e_targ.Data).';
    sim_output.attitude.targ_phi = rad2deg(e_targ.Data(:,1)); % deg
    sim_output.attitude.targ_theta = rad2deg(e_targ.Data(:,2)); % deg
    sim_output.attitude.targ_psi = rad2deg(e_targ.Data(:,3)); % deg
    
    sim_output.attitude.trueerr_phi = rad2deg(e_trueerr.Data(:,1)); % deg
    sim_output.attitude.trueerr_theta = rad2deg(e_trueerr.Data(:,2)); % deg
    sim_output.attitude.trueerr_psi = rad2deg(e_trueerr.Data(:,3)); % deg
    
    sim_output.attitude.A = permute(A.Data, [3,1,2]); % unitless, make time the first index
    sim_output.attitude.princ2inert = permute(A.Data, [3,2,1]); % inverse of A, rotates principal to inertial
    
    sim_output.attitude.A_est =  permute(A_est.Data, [3,1,2]);
    sim_output.attitude.A_esterr =  permute(A_esterr.Data, [3,1,2]);
    sim_output.attitude.A_trueerr = permute(A_trueerr.Data, [3,1,2]);
    sim_output.attitude.A_err = permute(A_err.Data, [3,1,2]); % unitless, make time the first index, rotates actual to target
    sim_output.attitude.A_targ = permute(A_targ.Data, [3,1,2]); % unitless, make time the first index, rotates inertial to target
    %% For attitude plots in inertial axes
    sim_output.attitude.w1 = zeros(size(sim_output.attitude.wx));
    sim_output.attitude.w2 = zeros(size(sim_output.attitude.wx));
    sim_output.attitude.w3 = zeros(size(sim_output.attitude.wx));
    
    sim_output.attitude.L1 = zeros(size(sim_output.attitude.Lx));
    sim_output.attitude.L2 = zeros(size(sim_output.attitude.Lx));
    sim_output.attitude.L3 = zeros(size(sim_output.attitude.Lx));
    
    W = [sim_output.attitude.wx, sim_output.attitude.wy, sim_output.attitude.wz];
    
    for i = 1:size(sim_output.attitude.A,1)
       R = squeeze(sim_output.attitude.princ2inert(i,:,:)); 
       Winertial = (R*(W(i,:).')).'; % convert to inertial axes
       sim_output.attitude.w1(i) = Winertial(1);
       sim_output.attitude.w2(i) = Winertial(2);
       sim_output.attitude.w3(i) = Winertial(3);
       L = (R*(diag(sim_constants.I_princ)*(W(i,:).') ... % calculate momentum in principal, convert to inertial
                + sim_output.attitude.L_r(i)*sim_constants.r_rotor)).'; 
       sim_output.attitude.L1(i) = L(1);
       sim_output.attitude.L2(i) = L(2);
       sim_output.attitude.L3(i) = L(3);
    end                        
end
