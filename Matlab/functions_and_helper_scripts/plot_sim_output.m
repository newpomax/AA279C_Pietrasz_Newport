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

    plot(downsampled_time, downsample(sim_output.attitude.wx, df));

    title_text = ['w_x of ', mission_name, ' in s/c principle axes'];
    title(title_text);
    xlabel(time_label);
    ytickformat('%.3g');
    ylabel('w_x [rad/s]');
    
    % w_y
    subplot(2,3,2);
    hold on;

    plot(downsampled_time, downsample(sim_output.attitude.wy, df));

    title_text = ['w_y of ', mission_name, ' in s/c principle axes'];
    title(title_text);
    xlabel(time_label);
    ytickformat('%.3g');
    ylabel('w_y [rad/s]');
    
    % w_z
    subplot(2,3,3);
    hold on;

    plot(downsampled_time, downsample(sim_output.attitude.wz, df));

    title_text = ['w_z of ', mission_name, ' in s/c principle axes'];
    title(title_text);
    xlabel(time_label);
    ytickformat('%.3g');
    ylabel('w_z [rad/s]');
    
    % L_x
    subplot(2,3,4);
    hold on;

    plot(downsampled_time, downsample(sim_output.attitude.Lx, df));

    title_text = ['L_x of ', mission_name, ' in s/c principle axes'];
    title(title_text);
    xlabel(time_label);
    ytickformat('%.3g');
    ylabel('L_x [kg m^2 /s]');
    
    % L_y
    subplot(2,3,5);
    hold on;

    plot(downsampled_time, downsample(sim_output.attitude.wy, df));

    title_text = ['L_y of ', mission_name, ' in s/c principle axes'];
    title(title_text);
    xlabel(time_label);
    ytickformat('%.3g');
    ylabel('L_y [kg m^2 /s]');
    
    % L_z
    subplot(2,3,6);
    hold on;

    plot(downsampled_time, downsample(sim_output.attitude.wz, df));

    title_text = ['L_z of ', mission_name, ' in s/c principle axes'];
    title(title_text);
    xlabel(time_label);
    ytickformat('%.3g');
    ylabel('L_z [kg m^2 /s]');
    
    %% Plot s/c Torques in Principal Axes
    figure('Name', strcat(mission_name, ' Environmental Torques'));
    % M_x
    subplot(1,3,1);
    hold on;

    plot(downsampled_time, downsample(sim_output.ext_torques.Mx, df));

    title_text = ['M_x of ', mission_name, ' in s/c principle axes'];
    title(title_text);
    xlabel(time_label);
    ytickformat('%.3g');
    ylabel('M_x [Nm]');
    
    % M_y
    subplot(1,3,2);
    hold on;

    plot(downsampled_time, downsample(sim_output.ext_torques.My, df));

    title_text = ['M_y of ', mission_name, ' in s/c principle axes'];
    title(title_text);
    xlabel(time_label);
    ytickformat('%.3g');
    ylabel('M_y [Nm]');
    
    % M_z
    subplot(1,3,3);
    hold on;

    plot(downsampled_time, downsample(sim_output.ext_torques.Mz, df));

    title_text = ['M_z of ', mission_name, ' in s/c principle axes'];
    title(title_text);
    xlabel(time_label);
    ytickformat('%.3g');
    ylabel('M_z [Nm]');
    
    if plot_format.ext_torques_verbose
        plot_ext_torques(sim_constants, sim_output, plot_format);
    else
        figure('Name','Torque Comparison'); 
        
        semilogy(downsampled_time, downsample(vecnorm( ...
             [sim_output.ext_torques.grav_Mx, sim_output.ext_torques.grav_My, sim_output.ext_torques.grav_Mz].'...
            ), df), 'DisplayName','Gravitational');
        hold on;
        semilogy(downsampled_time, downsample(vecnorm( ...
            [sim_output.ext_torques.drag_Mx, sim_output.ext_torques.drag_My, sim_output.ext_torques.drag_Mz].'...
            ), df), 'DisplayName','Drag');
        semilogy(downsampled_time, downsample(vecnorm( ...
            [sim_output.ext_torques.mag_Mx, sim_output.ext_torques.mag_My, sim_output.ext_torques.mag_Mz].'...
            ), df), 'DisplayName','Magnetic Torque');
        semilogy(downsampled_time, downsample(vecnorm( ...
            [sim_output.ext_torques.SRP_Mx, sim_output.ext_torques.SRP_My, sim_output.ext_torques.SRP_Mz].'...
            ), df), 'DisplayName','SRP');
        
        xlabel(time_label); ylabel('Torque, N\cdot m');
        title_text = ['|M| of ', plot_format.mission_name, ' in s/c principle axes'];
        title(title_text); legend('location','best');
    end
    
    %% Plot Angular Velocity + Momentum in Inertial Axes
    figure('Name',strcat(mission_name, ' Inertial Axes Ang. Vel. & Mom.'));
    % w_1
    subplot(2,3,1);
    hold on;

    plot(downsampled_time, downsample(sim_output.attitude.w1, df));

    title_text = ['\omega_1 of ', mission_name, ' in inertial axes'];
    title(title_text);
    xlabel(time_label);
    ytickformat('%.3g');
    ylabel('\omega_1 [rad/s]');
    
    % w_2
    subplot(2,3,2);
    hold on;

    plot(downsampled_time, downsample(sim_output.attitude.w2, df));

    title_text = ['\omega_2 of ', mission_name, ' in inertial axes'];
    title(title_text);
    xlabel(time_label);
    ytickformat('%.3g');
    ylabel('\omega_2 [rad/s]');
    
    % w_3
    subplot(2,3,3);
    hold on;

    plot(downsampled_time, downsample(sim_output.attitude.w3, df));

    title_text = ['\omega_3 of ', mission_name, ' in inertial axes'];
    title(title_text);
    xlabel(time_label);
    ytickformat('%.3g');
    ylabel('\omega_3 [rad/s]');
    
    % L_1
    subplot(2,3,4);
    hold on;

    plot(downsampled_time, downsample(sim_output.attitude.L1, df));

    title_text = ['L_1 of ', mission_name, ' in inertial axes'];
    title(title_text);
    xlabel(time_label);
    ytickformat('%.3g');
    ylabel('L_1 [kg m^2 /s]');
    
    % L_2
    subplot(2,3,5);
    hold on;

    plot(downsampled_time, downsample(sim_output.attitude.L2, df));

    title_text = ['L_2 of ', mission_name, ' in inertial axes'];
    title(title_text);
    xlabel(time_label);
    ytickformat('%.3g');
    ylabel('L_2 [kg m^2 /s]');
    
    % L_3
    subplot(2,3,6);
    hold on;

    plot(downsampled_time, downsample(sim_output.attitude.L3, df));

    title_text = ['L_3 of ', mission_name, ' in inertial axes'];
    title(title_text);
    xlabel(time_label);
    ytickformat('%.3g');
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
    ytickformat('%.3g');
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
    ytickformat('%.3g');
    ylabel('\omega_R [rad/s]');
    
    subplot(3,2,4); hold on;
    scatter(downsampled_time, downsample(W(:,2),df), 2);
    title_text = ['\omega_T of ', mission_name];
    title(title_text);
    xlabel(time_label);
    ytickformat('%.3g');
    ylabel('\omega_T [rad/s]');
    
    subplot(3,2,6); hold on;
    scatter(downsampled_time, downsample(W(:,3),df), 2);
    title_text = ['\omega_N of ', mission_name];
    title(title_text);
    xlabel(time_label);
    ytickformat('%.3g');
    ylabel('\omega_N [rad/s]');
    
    %% Plot coordinate triads (for first orbit only)
    if plot_format.plot_triads
        plot_coordtriads(sim_constants, sim_output, plot_format);
    end

    %% Plot 312 Euler angles
    
    % Force phi between -180 and 180:
    phi = sim_output.attitude.phi;
    idx = phi > 180;
    phi(idx) = phi(idx) - 360;
    
    figure('Name',strcat(mission_name, ' 312 Euler Angles')); 
    hold on;
    subplot(1,3,1); hold on;
    plot(sim_output.time, phi);
    ytickformat('%.3g');
    ylabel('\phi [deg]');
    xlabel(time_label);
    
    % Force theta between -90 and 90:
    % Assumes no values btwn 90 and 270 (which there shouldn't be).
    theta = sim_output.attitude.theta;
    idx = theta > 270;
    theta(idx) = theta(idx) - 360;
    
    subplot(1,3,2); hold on;
    plot(sim_output.time, theta);
    ytickformat('%.3g');
    ylabel('\theta [deg]');
    xlabel(time_label);
    
    % Force phi between -180 and 180:
    psi = sim_output.attitude.psi;
    idx = psi > 180;
    psi(idx) = psi(idx) - 360;
    
    subplot(1,3,3); hold on;
    plot(sim_output.time, psi);
    ytickformat('%.3g');
    ylabel('\psi [deg]');
    xlabel(time_label);
    
    title_text = ['312 Euler Angles of ', mission_name];
    sgtitle(title_text);
    
    %% Plot attitude control error
    
    % Force phi between -180 and 180:
    phi = sim_output.attitude.err_phi;
    idx = phi > 180;
    phi(idx) = phi(idx) - 360;
    phi_targ = sim_output.attitude.targ_phi;
    idx = phi_targ > 180;
    phi_targ(idx) = phi_targ(idx) - 360;
    
    figure('Name',strcat(mission_name, ' Attitude Error')); 
    hold on;
    subplot(1,3,1); hold on;
    plot(sim_output.time, phi_targ);
    plot(sim_output.time, phi);
    ytickformat('%.3g');
    ylabel('\phi [deg]');
    xlabel(time_label);
    legend('Target','Error');
    
    % Force theta between -90 and 90:
    % Assumes no values btwn 90 and 270 (which there shouldn't be).
    theta = sim_output.attitude.err_theta;
    idx = theta > 270;
    theta(idx) = theta(idx) - 360;
    theta_targ = sim_output.attitude.targ_theta;
    idx = theta_targ > 270;
    theta_targ(idx) = theta_targ(idx) - 360;
    
    subplot(1,3,2); hold on;
    plot(sim_output.time, theta_targ);
    plot(sim_output.time, theta);
    ytickformat('%.3g');
    ylabel('\theta [deg]');
    xlabel(time_label);
    legend('Target','Error');
    
    % Force phi between -180 and 180:
    psi = sim_output.attitude.err_psi;
    idx = psi > 180;
    psi(idx) = psi(idx) - 360;
    psi_targ = sim_output.attitude.targ_psi;
    idx = psi_targ > 180;
    psi_targ(idx) = psi_targ(idx) - 360;
    
    subplot(1,3,3); hold on;
    plot(sim_output.time, psi_targ);
    plot(sim_output.time, psi);
    ytickformat('%.3g');
    ylabel('\psi [deg]');
    xlabel(time_label);
    legend('Target', 'Error');
    
    title_text = ['312 Euler Angles for Attitude Control Error ', mission_name];
    sgtitle(title_text);
    
end