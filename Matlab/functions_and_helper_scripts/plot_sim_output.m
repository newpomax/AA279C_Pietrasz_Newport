function plot_sim_output(sim_constants, sim_output, plot_format)
    % Plots simulation outputs using data extracted as in extract_sim_output.
    % Generates:
    %   - orbital elements over time
    %   - orbital element rates over time
    %   - optionally, a visualization of the orbit propagation. If you are 
    % simulating many orbits, this last plot becomes too crowded 
    %   - angular velocity and momentum over time
    %       - In s/c principal axes
    %       - In ECI axes
    %   - attitude herpolhode, RTN frame over time
    %   - attitude triads over time
    %   - Euler angles (312)

    mission_name = plot_format.mission_name;
    time_label = ['Time [', num2str(plot_format.time_increments), ']'];

    df = plot_format.downsample_factor;
    downsampled_time = downsample(sim_output.time, df);
    
    if plot_format.dock_plots
       set(0,'DefaultFigureWindowStyle','docked'); 
    end

    plot_orbital_elements(sim_constants, sim_output, plot_format);
    
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
    
    %% Plot s/c Torques in Principal Axes
    figure('Name',strcat(mission_name, ' Environmental Torques'));
    % w_x
    subplot(1,3,1);
    hold on;

    scatter(downsampled_time, downsample(sim_output.ext_torques.Mx, df), 2);

    title_text = ['M_x of ', mission_name, ' in s/c principle axes'];
    title(title_text);
    xlabel(time_label);
    ylabel('w_x [rad/s]');
    
    % w_y
    subplot(1,3,2);
    hold on;

    scatter(downsampled_time, downsample(sim_output.ext_torques.My, df), 2);

    title_text = ['w_y of ', mission_name, ' in s/c principle axes'];
    title(title_text);
    xlabel(time_label);
    ylabel('w_y [rad/s]');
    
    % w_z
    subplot(1,3,3);
    hold on;

    scatter(downsampled_time, downsample(sim_output.ext_torques.Mz, df), 2);

    title_text = ['w_z of ', mission_name, ' in s/c principle axes'];
    title(title_text);
    xlabel(time_label);
    ylabel('w_z [rad/s]');
    
    % Maybe control torques here?
    
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
    
    %% Plot Herpolhode & RTN Frame Ang Vel
    figure('Name',strcat(mission_name, ' Herpolhode & RTN'));
    subplot(3,2,[1 3 5]);
    hold on; 
    
    plot3(sim_output.attitude.w1,sim_output.attitude.w2,sim_output.attitude.w3)
    L = 100*norm([sim_output.attitude.L1(1) sim_output.attitude.L2(1) sim_output.attitude.L3(1)]);
    
    start = mean([sim_output.attitude.w1 sim_output.attitude.w2 sim_output.attitude.w3]);
    finish = [sim_output.attitude.L1(1) sim_output.attitude.L2(1) sim_output.attitude.L3(1)]/L;
    quiver3(start(1),start(2),start(3),finish(1),finish(2),finish(3),'LineWidth',2);
  
    title_text = ['Herpohlode of ', mission_name];
    title(title_text);
    xlabel('\omega_1');
    ylabel('\omega_2');
    zlabel('\omega_3');
    
    legend('Herpolhode','Ang. Momentum');
    view(3);
    axis equal;
    pbaspect([1 1 1]);
    
    % Transform omega vector from ECI to RTN
    W = zeros(length(sim_output.attitude.w1),3);
    for i = 1:size(sim_output.positions.RTN2ECI,1)
        W(i,:) = ((squeeze(sim_output.positions.RTN2ECI(i,:,:)).')*([sim_output.attitude.w1(i);sim_output.attitude.w2(i);sim_output.attitude.w3(i)])).';
    end
    subplot(3,2,2); hold on;
    scatter(downsampled_time, downsample(W(:,1),df), 2);
    title_text = ['\omega_R of ', mission_name];
    title(title_text);
    xlabel(time_label);
    ylabel('\omega_R [rad/s]');
    
    subplot(3,2,4); hold on;
    scatter(downsampled_time, downsample(W(:,2),df), 2);
    title_text = ['\omega_T of ', mission_name];
    title(title_text);
    xlabel(time_label);
    ylabel('\omega_T [rad/s]');
    
    subplot(3,2,6); hold on;
    scatter(downsampled_time, downsample(W(:,3),df), 2);
    title_text = ['\omega_N of ', mission_name];
    title(title_text);
    xlabel(time_label);
    ylabel('\omega_N [rad/s]');
    
    %% Plot coordinate triads (for first orbit only)
    if plot_format.plot_triads
        plot_coordtriads(sim_constants, sim_output, plot_format);
    end

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
    
end