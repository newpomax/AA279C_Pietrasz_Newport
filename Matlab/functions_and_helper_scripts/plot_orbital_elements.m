function plot_orbital_elements(sim_constants, sim_output, plot_format)
    % Subfunction of plot_sim_output that plots all the orbital elements
    % and an orbital visual if desired.

    mission_name = plot_format.mission_name;
    time_label = ['Time [', num2str(plot_format.time_increments), ']'];

    df = plot_format.downsample_factor;
    downsampled_time = downsample(sim_output.time, df);
    
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
        plot3(sim_output.positions.x_XYZ, sim_output.positions.y_XYZ, ...
            sim_output.positions.z_XYZ, 'LineWidth', 1);

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

end