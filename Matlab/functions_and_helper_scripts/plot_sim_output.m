function plot_sim_output(sim_constants, sim_output, plot_format)
    % Plots simulation outputs using data extracted as in extract_sim_output.
    % Generates 2 figures, plus one optional figure (depending on whether
    % plot_format.plot_orbit_visuals is true or false): 
    %   - orbital elements over time
    %   - orbital element rates over time
    %   - optionally, a visualization of the orbit propagation. If you are 
    % simulating many orbits, this last plot becomes too crowded 

    mission_name = plot_format.mission_name;
    time_label = ['Time [', num2str(plot_format.time_increments), ']'];

    df = plot_format.downsample_factor;
    downsampled_time = downsample(sim_output.time, df);
    
    if plot_format.dock_plots
       set(0,'DefaultFigureWindowStyle','docked'); 
    end
    
    if plot_format.plot_orbit_visuals
        figure('Name', strcat(mission_name, ' orbit visualization'));
        if ~plot_format.dock_plots
            set(gcf, 'Position',  [0, 0, 500, 800]);
            subplot(2,1,1);
        else
            subplot(1,3,1);
        end

        % 3D
        hold on;

        % Plot the Earth for context
        R_Earth = sim_constants.R_Earth;
        [xE, yE, zE] = ellipsoid(0, 0, 0, R_Earth, R_Earth, R_Earth, 20);
        surface(xE, yE, zE, 'FaceColor', 'blue', 'EdgeColor', 'black'); 

        % Plot Simulink output
        plot3(sim_output.positions.x_XYZ, sim_output.positions.y_XYZ, sim_output.positions.z_XYZ, ...
            'LineWidth', 1);

        % label n things
        title_text = ['Orbit of ', mission_name, ...
            ' in ECEF reference frame'];
        title(title_text)
        xlabel('km');
        ylabel('km');
        zlabel('km');
        axis equal;
        pbaspect([1 1 1]);
        view(3);
        grid on;

        % Groundtrack
        if plot_format.dock_plots
            subplot(1,3,[2 3]);
        else
            subplot(2,1,2);
        end
        hold on;

        % Plot the Earth for context
        geoshow('landareas.shp', 'FaceColor', [0.5 1.0 0.5]);

        % Plot Simulink output
        geoshow(sim_output.positions.latitude_d, sim_output.positions.longitude_d);

        % label n things
        title_text = ['Groundtrack of ', plot_format.mission_name, ' ', ...
            ' in the Mercator lat/long projection'];
        title(title_text);
        xlabel('lat');
        ylabel('long');
        grid on;
    end

    % OEs
    figure('Name', strcat(mission_name, ' OEs'));
    if ~plot_format.dock_plots
        set(gcf, 'Position',  [0, 0, 1200, 600])
    end

    % semimajor axis
    subplot(2,3,1);
    hold on;

    scatter(downsampled_time, downsample(sim_output.OE.a, df), 2);

    title_text = ['Semimajor axis of ', mission_name];
    title(title_text);
    xlabel(time_label);
    ylabel('a [km]');

    % eccentricity
    subplot(2,3,4);
    hold on;

    scatter(downsampled_time, downsample(sim_output.OE.e, df), 2);

    title_text = ['Eccentricity of ', mission_name];
    title(title_text);
    xlabel(time_label);
    ylabel('e');

    % inclination
    subplot(2,3,3);
    hold on;

    scatter(downsampled_time, downsample(sim_output.OE.i, df), 2);

    title_text = ['Inclination of ', mission_name];
    title(title_text);
    xlabel(time_label);
    ylabel('i [deg]');

    % RAAN
    subplot(2,3,2);
    hold on;

    scatter(downsampled_time, downsample(sim_output.OE.RAAN, df), 2);

    title_text = ['RAAN of ', mission_name];
    title(title_text);
    xlabel(time_label);
    ylabel('RAAN [deg]');

    % w
    subplot(2,3,5);
    hold on;

    scatter(downsampled_time, downsample(sim_output.OE.w, df), 2);

    title_text = ['Argument of Periapsis of ', mission_name];
    title(title_text);
    xlabel(time_label);
    ylabel('w [deg]');

    % E
    subplot(2,3,6);
    hold on;

    scatter(downsampled_time, downsample(sim_output.OE.E, df), 2);

    title_text = ['Eccentric anomaly of ', mission_name];
    title(title_text);
    xlabel(time_label);
    ylabel('E [rad]');

    %% Plot dOE/dt
    
    figure('Name', strcat(mission_name, ' dOE/dts'));
    if ~plot_format.dock_plots
        set(gcf, 'Position',  [1200, 0, 1200, 600])
    end

    % semimajor axis
    subplot(2,3,1);
    hold on;

    scatter(downsampled_time, ...
        downsample(sim_output.dOE.a, df).*1000, 2);

    title_text = ['da/dt of ', mission_name];
    title(title_text);
    xlabel(time_label);
    ylabel('da/dt [m/s]');

    % de/dt
    subplot(2,3,4);
    hold on;

    scatter(downsampled_time, downsample(sim_output.dOE.e, df), 2);

    title_text = ['de/dt of ', mission_name];
    title(title_text);
    xlabel(time_label);
    ylabel('de/dt');

    % di/dt
    subplot(2,3,3);
    hold on;

    scatter(downsampled_time, ...
        downsample(sim_output.dOE.i, df)*10^3, 2);

    title_text = ['di/dt of ', mission_name];
    title(title_text);
    xlabel(time_label);
    ylabel('di/dt [mdeg/sec]');

    % dRAAN/dt
    subplot(2,3,2);
    hold on;
    scatter(downsampled_time, ...
        downsample(sim_output.dOE.RAAN, df).*10^3, 2); 

    title_text = ['dRAAN/dt of ', mission_name];
    title(title_text);
    xlabel(time_label);
    ylabel('dRAAN/dt [mdeg/s]');

    % dw/dt
    subplot(2,3,5);
    hold on;

    scatter(downsampled_time, downsample(sim_output.dOE.w, df), 2);

    title_text = ['dw/dt of ', mission_name];
    title(title_text);
    xlabel(time_label);
    ylabel('dw/dt [deg/s]');

    % dE/dt
    subplot(2,3,6);
    hold on;

    scatter(downsampled_time, downsample(sim_output.dOE.E, df), 2);

    title_text = ['dE/dt of ', mission_name, ' relative to n'];
    title(title_text);
    xlabel(time_label);
    ylabel('dE/dt [rad/s]');
    
    %% Plot Angular Velocity + Momentum in Principal Axes
    figure('Name',strcat(mission_name, ' Princ. Axes Ang. Vel. & Mom.'));
    % w_x
    subplot(2,3,1);
    hold on;

    scatter(downsampled_time, downsample(sim_output.attitude.wx, df), 2);

    title_text = ['w_x of ', mission_name, ' in s/c principle axes'];
    title(title_text);
    xlabel(time_label);
    ylabel('w_x [rad/s]');
    
    % w_y
    subplot(2,3,2);
    hold on;

    scatter(downsampled_time, downsample(sim_output.attitude.wy, df), 2);

    title_text = ['w_y of ', mission_name, ' in s/c principle axes'];
    title(title_text);
    xlabel(time_label);
    ylabel('w_y [rad/s]');
    
    % w_z
    subplot(2,3,3);
    hold on;

    scatter(downsampled_time, downsample(sim_output.attitude.wz, df), 2);

    title_text = ['w_z of ', mission_name, ' in s/c principle axes'];
    title(title_text);
    xlabel(time_label);
    ylabel('w_z [rad/s]');
    
    % L_x
    subplot(2,3,4);
    hold on;

    scatter(downsampled_time, downsample(sim_output.attitude.Lx, df), 2);

    title_text = ['L_x of ', mission_name, ' in s/c principle axes'];
    title(title_text);
    xlabel(time_label);
    ylabel('L_x [kg m^2 /s]');
    
    % L_y
    subplot(2,3,5);
    hold on;

    scatter(downsampled_time, downsample(sim_output.attitude.wy, df), 2);

    title_text = ['L_y of ', mission_name, ' in s/c principle axes'];
    title(title_text);
    xlabel(time_label);
    ylabel('L_y [kg m^2 /s]');
    
    % L_z
    subplot(2,3,6);
    hold on;

    scatter(downsampled_time, downsample(sim_output.attitude.wz, df), 2);

    title_text = ['L_z of ', mission_name, ' in s/c principle axes'];
    title(title_text);
    xlabel(time_label);
    ylabel('L_z [kg m^2 /s]');
    
     %% Plot Angular Velocity + Momentum in Inertial Axes
    figure('Name',strcat(mission_name, ' Inertial Axes Ang. Vel. & Mom.'));
    % w_1
    subplot(2,3,1);
    hold on;

    scatter(downsampled_time, downsample(sim_output.attitude.w1, df), 2);

    title_text = ['\omega_1 of ', mission_name, ' in inertial axes'];
    title(title_text);
    xlabel(time_label);
    ylabel('\omega_1 [rad/s]');
    
    % w_2
    subplot(2,3,2);
    hold on;

    scatter(downsampled_time, downsample(sim_output.attitude.w2, df), 2);

    title_text = ['\omega_2 of ', mission_name, ' in inertial axes'];
    title(title_text);
    xlabel(time_label);
    ylabel('\omega_2 [rad/s]');
    
    % w_3
    subplot(2,3,3);
    hold on;

    scatter(downsampled_time, downsample(sim_output.attitude.w3, df), 2);

    title_text = ['\omega_3 of ', mission_name, ' in inertial axes'];
    title(title_text);
    xlabel(time_label);
    ylabel('\omega_3 [rad/s]');
    
    % L_1
    subplot(2,3,4);
    hold on;

    scatter(downsampled_time, downsample(sim_output.attitude.L1, df), 2);

    title_text = ['L_1 of ', mission_name, ' in inertial axes'];
    title(title_text);
    xlabel(time_label);
    ylabel('L_1 [kg m^2 /s]');
    
    % L_2
    subplot(2,3,5);
    hold on;

    scatter(downsampled_time, downsample(sim_output.attitude.L2, df), 2);

    title_text = ['L_2 of ', mission_name, ' in inertial axes'];
    title(title_text);
    xlabel(time_label);
    ylabel('L_2 [kg m^2 /s]');
    
    % L_3
    subplot(2,3,6);
    hold on;

    scatter(downsampled_time, downsample(sim_output.attitude.L3, df), 2);

    title_text = ['L_3 of ', mission_name, ' in inertial axes'];
    title(title_text);
    xlabel(time_label);
    ylabel('L_3 [kg m^2 /s]');
    
    %% Plot Herpohlode & RTN Frame Ang Vel
    figure('Name',strcat(mission_name, ' Herpolhode & RTN'));
    subplot(1,2,1);
    hold on; 
    
    plot3(sim_output.attitude.w1,sim_output.attitude.w2,sim_output.attitude.w3)
    L = 100*norm([sim_output.attitude.L1(1) sim_output.attitude.L2(1) sim_output.attitude.L3(1)]);
    
    start = mean([sim_output.attitude.w1 sim_output.attitude.w2 sim_output.attitude.w3]);
    finish = [sim_output.attitude.L1(1) sim_output.attitude.L2(1) sim_output.attitude.L3(1)]/L;
    quiver3(start(1),start(2),start(3),finish(1),finish(2),finish(3),'LineWidth',2);

    legend('Herpolhode','Ang. Momentum');
    view(3);
    axis equal;
    
    title_text = ['Herpohlode of ', mission_name];
    title(title_text);
    xlabel('\omega_1');
    ylabel('\omega_2');
    zlabel('\omega_3');
    
    subplot(1,2,2); hold on;
    % Transform omega vector from ECI to RTN
    W = zeros(length(sim_output.attitude.w1),3);
    for i = 1:size(sim_output.positions.RTN2ECI,1)
        W(i,:) = ((squeeze(sim_output.positions.RTN2ECI(i,:,:)).')*([sim_output.attitude.w1(i);sim_output.attitude.w2(i);sim_output.attitude.w3(i)])).';
    end
    plot3(W(:,1),W(:,2),W(:,3));
    view(3);
    axis equal;
    
    title_text = ['RTN Frame Ang Velocity of ', mission_name];
    title(title_text);
    xlabel('\omega_R');
    ylabel('\omega_T');
    zlabel('\omega_N');
    
    %% Plot 312 Euler angles
    figure('Name',strcat(mission_name, ' 312 Euler Angles')); 
    hold on;
    subplot(1,3,1); hold on;
    plot(sim_output.time,  sim_output.attitude.phi);
    ylabel('\phi [deg]');
    xlabel(time_label);
    
    subplot(1,3,2); hold on;
    plot(sim_output.time,  sim_output.attitude.theta);
    ylabel('\theta [deg]');
    xlabel(time_label);
    
    subplot(1,3,3); hold on;
    plot(sim_output.time,  sim_output.attitude.psi);
    ylabel('\psi [deg]');
    xlabel(time_label);
    
    title_text = ['312 Euler Angles of ', mission_name];
    sgtitle(title_text);
    
    %% Plot coordinate triads (for first orbit only)
    if plot_format.plot_triads
        plot_coordtriads(sim_constants, sim_output, plot_format);
    end
end