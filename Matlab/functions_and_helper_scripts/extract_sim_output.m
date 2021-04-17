function sim_output = extract_sim_output(sim_constants, plot_format, ...
    OE, dOE_dt, w, q, A, ECEF_positions, geod_positions)
    %% Post-process data output from simulation for plotting and calcs

    sim_output = struct;
    sim_output.time = OE.time/plot_format.seconds_to_increment;

    %% Postions for orbit visualizations

    sim_output.positions.x_XYZ = ECEF_positions.Data(:,1);
    sim_output.positions.y_XYZ = ECEF_positions.Data(:,2);
    sim_output.positions.z_XYZ = ECEF_positions.Data(:,3);

    sim_output.positions.latitude_d = geod_positions.Data(:,1);
    sim_output.positions.longitude_d = geod_positions.Data(:,2);

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
    
    %% For attitude plots in principal axes
    
    sim_output.attitude.wx = w.Data(:,1); % rad s-1, principal axes
    sim_output.attitude.wy = w.Data(:,2); % rad s-1, principal axes
    sim_output.attitude.wz = w.Data(:,3); % rad s-1, principal axes
    
    sim_output.attitude.Lx = sim_constants.I_princ(1)*sim_output.attitude.wx; % kg m^2/s, principal axes
    sim_output.attitude.Ly = sim_constants.I_princ(2)*sim_output.attitude.wy; % kg m^2/s, principal axes
    sim_output.attitude.Lz = sim_constants.I_princ(3)*sim_output.attitude.wz; % kg m^2/s, principal axes
    
    sim_output.attitude.q1 = q.Data(:,1); % unitless
    sim_output.attitude.q2 = q.Data(:,2); % unitless
    sim_output.attitude.q3 = q.Data(:,3); % unitless
    sim_output.attitude.q4 = q.Data(:,4); % unitless
    
    sim_output.attitude.A = permute(A.Data, [3,1,2]); % unitless, make time the first index
    
    %% For attitude plots in inertial axes
    sim_output.attitude.princ2inert = zeros(size(sim_output.attitude.A));
    sim_output.attitude.w1 = zeros(size(sim_output.attitude.wx));
    sim_output.attitude.w2 = zeros(size(sim_output.attitude.wx));
    sim_output.attitude.w3 = zeros(size(sim_output.attitude.wx));
    sim_output.attitude.L1 = zeros(size(sim_output.attitude.Lx));
    sim_output.attitude.L2 = zeros(size(sim_output.attitude.Lx));
    sim_output.attitude.L3 = zeros(size(sim_output.attitude.Lx));
    
    W = [sim_output.attitude.wx, sim_output.attitude.wy, sim_output.attitude.wz];
    for i = 1:size(sim_output.attitude.A,1)
       R = squeeze(sim_output.attitude.A(i,:,:)).'; % inverse, rotates principal to inertial
       sim_output.attitude.princ2inert(i,:,:) = R;
       Winertial = (R*(W.')).'; % convert to inertial axes
       sim_output.attitude.w1 = Winertial(:,1);
       sim_output.attitude.w2 = Winertial(:,2);
       sim_output.attitude.w3 = Winertial(:,3);
       L = (R*(diag(sim_constants.I_princ)*(W.'))).'; % calculate momentum in principal, convert to inertial
       sim_output.attitude.L1 = L(:,1);
       sim_output.attitude.L2 = L(:,2);
       sim_output.attitude.L3 = L(:,3);
    end                                               
end
