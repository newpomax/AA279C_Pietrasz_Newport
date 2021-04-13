function sim_output = extract_sim_output(sim_constants, plot_format, ...
    OE, dOE_dt, w, ECEF_positions, geod_positions)
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
    sim_output.OE.w = OE.Data(:,5);; % deg
    sim_output.OE.E = OE.Data(:,6); % rad

    sim_output.dOE.a = dOE_dt.Data(:,1); % km
    sim_output.dOE.e = dOE_dt.Data(:,2); % deg s-1
    sim_output.dOE.i = dOE_dt.Data(:,3); % deg s-1

    sim_output.dOE.RAAN = dOE_dt.Data(:,4); % deg s-1
    sim_output.dOE.w = dOE_dt.Data(:,5); % deg s-1
    sim_output.dOE.E = dOE_dt.Data(:,6); % deg s-1
    
    % For attitude plots
    
    sim_output.w.x = w.Data(:,1); % km
    sim_output.w.y = w.Data(:,2); % deg
    sim_output.w.z = w.Data(:,3); % deg
    
end
