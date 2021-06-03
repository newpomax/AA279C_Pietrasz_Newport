function plot_batch_sims(sim_constants, sim_output, plot_format)

    mission_name = plot_format.mission_name;
    time_label = ['Time [', num2str(plot_format.time_increments), ']'];

    df = plot_format.downsample_factor;
    downsampled_time = downsample(sim_output.time, df);
    
    if plot_format.dock_plots
       set(0,'DefaultFigureWindowStyle','docked'); 
    end
    
%     %% Plot Angular Velocity + Momentum in Principal Axes
%     figure('Name',strcat(mission_name, ' Princ. Axes Ang. Vel. & Mom.'));
%     % w_x
%     subplot(2,3,1);
%     hold on;
% 
%     plot(downsampled_time, downsample(sim_output.attitude.wx, df));
% 
%     title_text = ['w_x of ', mission_name, ' in s/c principle axes'];
%     title(title_text);
%     xlabel(time_label);
%     ytickformat('%.3g');
%     ylabel('w_x [rad/s]');
%     
%     % w_y
%     subplot(2,3,2);
%     hold on;
% 
%     plot(downsampled_time, downsample(sim_output.attitude.wy, df));
% 
%     title_text = ['w_y of ', mission_name, ' in s/c principle axes'];
%     title(title_text);
%     xlabel(time_label);
%     ytickformat('%.3g');
%     ylabel('w_y [rad/s]');
%     
%     % w_z
%     subplot(2,3,3);
%     hold on;
% 
%     plot(downsampled_time, downsample(sim_output.attitude.wz, df));
% 
%     title_text = ['w_z of ', mission_name, ' in s/c principle axes'];
%     title(title_text);
%     xlabel(time_label);
%     ytickformat('%.3g');
%     ylabel('w_z [rad/s]');
%     
%     % L_x
%     subplot(2,3,4);
%     hold on;
% 
%     plot(downsampled_time, downsample(sim_output.attitude.Lx, df));
% 
%     title_text = ['L_x of ', mission_name, ' in s/c principle axes'];
%     title(title_text);
%     xlabel(time_label);
%     ytickformat('%.3g');
%     ylabel('L_x [kg m^2 /s]');
%     
%     % L_y
%     subplot(2,3,5);
%     hold on;
% 
%     plot(downsampled_time, downsample(sim_output.attitude.wy, df));
% 
%     title_text = ['L_y of ', mission_name, ' in s/c principle axes'];
%     title(title_text);
%     xlabel(time_label);
%     ytickformat('%.3g');
%     ylabel('L_y [kg m^2 /s]');
%     
%     % L_z
%     subplot(2,3,6);
%     hold on;
% 
%     plot(downsampled_time, downsample(sim_output.attitude.wz, df));
% 
%     title_text = ['L_z of ', mission_name, ' in s/c principle axes'];
%     title(title_text);
%     xlabel(time_label);
%     ytickformat('%.3g');
%     ylabel('L_z [kg m^2 /s]');

    %% Plot attitude control torques, magnetorquer, and momentum wheel info
    figure('Name',strcat(mission_name, ' Control Torque')); hold on;
    labs = {'x','y','z'};
    for i = 1:3
        subplot(1,3,i); hold on;
        title_text = ['M_{control,' labs{i} '} of ' mission_name];
        title(title_text);
        xlabel(time_label);
        ylabel('Torque, Nm');
        plot(downsampled_time,downsample(sim_output.attitude.M_c(:,i),df));
        plot(downsampled_time,downsample(sim_output.attitude.M_req(:,i),df));
        legend('True Control','ADCS Torque Request');
    end
    
%     figure('Name',strcat(mission_name, ' Magnetorquers')); hold on;
%     labs = {'x','y','z'};
%     for i = 1:3
%         subplot(1,3,i); hold on;
%         title_text = ['Magnetorquer Dipole Along ' labs{i} ' for ' mission_name];
%         title(title_text);
%         xlabel(time_label);
%         ylabel(['m_' labs{i} ', Am^2']);
%         plot(downsampled_time,downsample(sim_output.attitude.m_magtor(:,i),df));
%     end
%     
%     figure('Name',strcat(mission_name, ' Mom. Wheel')); hold on;
%     subplot(1,2,1); hold on;
%     title_text = ['Momentum Wheel Speed for ' mission_name];
%     title(title_text);
%     xlabel(time_label);
%     ylabel('\omega_r , rad/s');
%     plot(downsampled_time,downsample(sim_output.attitude.w_r,df));
%     subplot(1,2,2); hold on;
%     title_text = ['Momentum Wheel Accel. for ' mission_name];
%     title(title_text);
%     xlabel(time_label);
%     ylabel('d\omega_r/dt, rad/s^2');
%     plot(downsampled_time,downsample(sim_output.attitude.dw_rdt,df));
    
    %% Plot 312 Euler angles
    
    if false
        % Force phi between -180 and 180:
        phi = sim_output.attitude.phi;
        idx = phi > 180;
        phi(idx) = phi(idx) - 360;
        phi_est = sim_output.attitude.est_phi;
        idx = phi_est > 180;
        phi_est(idx) = phi_est(idx) - 360;

        figure('Name',strcat(mission_name, ' 312 Euler Angles')); 
        hold on;
        subplot(1,3,1); hold on;
        plot(sim_output.time, phi);
        plot(sim_output.time, phi_est);
        ytickformat('%.3g');
        ylabel('\phi [deg]');
        xlabel(time_label);
        legend('True','Estimate');

        % Force theta between -90 and 90:
        % Assumes no values btwn 90 and 270 (which there shouldn't be).
        theta = sim_output.attitude.theta;
        idx = theta > 270;
        theta(idx) = theta(idx) - 360;
        theta_est = sim_output.attitude.est_theta;
        idx = theta_est > 180;
        theta_est(idx) = theta_est(idx) - 360;

        subplot(1,3,2); hold on;
        plot(sim_output.time, theta);
        plot(sim_output.time, theta_est);
        ytickformat('%.3g');
        ylabel('\theta [deg]');
        xlabel(time_label);
        legend('True','Estimate');

        % Force phi between -180 and 180:
        psi = sim_output.attitude.psi;
        idx = psi > 180;
        psi(idx) = psi(idx) - 360;
        psi_est = sim_output.attitude.est_psi;
        idx = psi_est > 180;
        psi_est(idx) = psi_est(idx) - 360;

        subplot(1,3,3); hold on;
        plot(sim_output.time, psi);
        plot(sim_output.time, psi_est);
        ytickformat('%.3g');
        ylabel('\psi [deg]');
        xlabel(time_label);
        legend('True','Estimate');

        title_text = ['312 Euler Angles of ', mission_name];
        sgtitle(title_text);
    end
    
    %% Plot attitude control error
    
    % Force phi between -180 and 180:
    phi = sim_output.attitude.err_phi;
    idx = phi > 180;
    phi(idx) = phi(idx) - 360;
    phi_targ = sim_output.attitude.targ_phi;
    idx = phi_targ > 180;
    phi_targ(idx) = phi_targ(idx) - 360;
    phi_true = sim_output.attitude.trueerr_phi;
    idx = phi_true > 180;
    phi_true(idx) = phi_true(idx) - 360;
    
    figure('Name', strcat(mission_name, 'Target Attitude Error')); 
    hold on;
    subplot(1,3,1); hold on;
    plot(downsampled_time, downsample(phi_targ,df));
    plot(downsampled_time, downsample(phi,df));
    plot(downsampled_time, downsample(phi_true,df));
    ytickformat('%.3g');
    ylabel('\phi [deg]');
    xlabel(time_label);
    legend('Target','Est. Error','True Error');
    
    % Force theta between -90 and 90:
    % Assumes no values btwn 90 and 270 (which there shouldn't be).
    theta = sim_output.attitude.err_theta;
    idx = theta > 270;
    theta(idx) = theta(idx) - 360;
    theta_targ = sim_output.attitude.targ_theta;
    idx = theta_targ > 270;
    theta_targ(idx) = theta_targ(idx) - 360;
    theta_true = sim_output.attitude.trueerr_theta;
    idx = theta_true > 180;
    theta_true(idx) = theta_true(idx) - 360;
    
    subplot(1,3,2); hold on;
    plot(downsampled_time, downsample(theta_targ,df));
    plot(downsampled_time, downsample(theta,df));
    plot(downsampled_time, downsample(theta_true,df));
    ytickformat('%.3g');
    ylabel('\theta [deg]');
    xlabel(time_label);
    legend('Target','Est. Error','True Error');
    
    % Force phi between -180 and 180:
    psi = sim_output.attitude.err_psi;
    idx = psi > 180;
    psi(idx) = psi(idx) - 360;
    psi_targ = sim_output.attitude.targ_psi;
    idx = psi_targ > 180;
    psi_targ(idx) = psi_targ(idx) - 360;
    psi_true = sim_output.attitude.trueerr_psi;
    idx = psi_true > 180;
    psi_true(idx) = psi_true(idx) - 360;
    
    subplot(1,3,3); hold on;
    plot(downsampled_time, downsample(psi_targ,df));
    plot(downsampled_time, downsample(psi,df));
    plot(downsampled_time, downsample(psi_true,df));
    ytickformat('%.3g');
    ylabel('\psi [deg]');
    xlabel(time_label);
    legend('Target', 'Est. Error', 'True Error');
    
    title_text = ['312 Euler Angles for Attitude Control Error ', mission_name];
    sgtitle(title_text);
    
    %% Plot Estimation Error
    
    if false
        % Force phi between -180 and 180:
        phi_esterr = sim_output.attitude.esterr_phi;
        idx = phi_esterr > 180;
        phi_esterr(idx) = phi_esterr(idx) - 360;

        figure('Name',strcat(mission_name, ' Attitude Estimate')); 
        hold on;
        subplot(1,3,1); hold on;
        plot(downsampled_time, downsample(phi_esterr,df));
        ytickformat('%.3g');
        ylabel('\phi [deg]');
        xlabel(time_label);

        % Force theta between -90 and 90:
        % Assumes no values btwn 90 and 270 (which there shouldn't be).
        theta_esterr = sim_output.attitude.esterr_theta;
        idx = theta_esterr > 180;
        theta_esterr(idx) = theta_esterr(idx) - 360;

        subplot(1,3,2); hold on;
        plot(downsampled_time, downsample(theta_esterr,df));
        ytickformat('%.3g');
        ylabel('\theta [deg]');
        xlabel(time_label);

        % Force phi between -180 and 180:
        psi_esterr = sim_output.attitude.esterr_psi;
        idx = psi_esterr > 180;
        psi_esterr(idx) = psi_esterr(idx) - 360;

        subplot(1,3,3); hold on;
        plot(downsampled_time, downsample(psi_esterr,df));
        ytickformat('%.3g');
        ylabel('\psi [deg]');
        xlabel(time_label);

        title_text = ['Attitude Estimation Error of ', mission_name];
        sgtitle(title_text);
    end
    
    %% Plot SC Mode
    figure('Name',strcat(mission_name, ' Attitude Mode')); hold on;
    title_text = ['Spacecraft Attitude Control Mode for ', mission_name];
    title(title_text);
    xlabel(time_label);
    plot(sim_output.time, sim_output.attitude.SC_Mode);
    ylabel('SC Mode');
    [~,y] = enumeration('SCModeEnum');
    yticklabels(y); yticks([0 1 2 3]);

end

