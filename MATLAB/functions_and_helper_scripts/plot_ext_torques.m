function plot_ext_torques(sim_constants, sim_output, plot_format)
    % Subfunction of plot_sim_output that plots all the external s/c
    % torques.
    
    mission_name = plot_format.mission_name;
    time_label = ['Time [', num2str(plot_format.time_increments), ']'];

    df = plot_format.downsample_factor;
    downsampled_time = downsample(sim_output.time, df);
    
    % Drag
    figure('Name',strcat(mission_name, ' Drag Torques'));
    % M_x
    subplot(1,3,1);
    hold on;

    plot(downsampled_time, downsample(sim_output.ext_torques.drag_Mx, df));

    title_text = ['Drag M_x of ', mission_name, ' in s/c principle axes'];
    title(title_text);
    xlabel(time_label);
    ytickformat('%.3g');
    ylabel('M_x [Nm]');

    % M_y
    subplot(1,3,2);
    hold on;

    plot(downsampled_time, downsample(sim_output.ext_torques.drag_My, df));

    title_text = ['Drag M_y of ', mission_name, ' in s/c principle axes'];
    title(title_text);
    xlabel(time_label);
    ytickformat('%.3g');
    ylabel('M_y [Nm]');

    % M_z
    subplot(1,3,3);
    hold on;

    plot(downsampled_time, downsample(sim_output.ext_torques.drag_Mz, df));

    title_text = ['Drag M_z of ', mission_name, ' in s/c principle axes'];
    title(title_text);
    xlabel(time_label);
    ytickformat('%.3g');
    ylabel('M_z [Nm]');

    % SRP
    figure('Name',strcat(mission_name, ' SRP Torques'));
    % M_x
    subplot(1,3,1);
    hold on;

    plot(downsampled_time, downsample(sim_output.ext_torques.SRP_Mx, df));

    title_text = ['SRP M_x of ', mission_name, ' in s/c principle axes'];
    title(title_text);
    xlabel(time_label);
    ytickformat('%.3g');
    ylabel('M_x [Nm]');

    % M_y
    subplot(1,3,2);
    hold on;

    plot(downsampled_time, downsample(sim_output.ext_torques.SRP_My, df));

    title_text = ['SRP M_y of ', mission_name, ' in s/c principle axes'];
    title(title_text);
    xlabel(time_label);
    ytickformat('%.3g');
    ylabel('M_y [Nm]');

    % M_z
    subplot(1,3,3);
    hold on;

    plot(downsampled_time, downsample(sim_output.ext_torques.SRP_Mz, df));

    title_text = ['SRP M_z of ', mission_name, ' in s/c principle axes'];
    title(title_text);
    xlabel(time_label);
    ytickformat('%.3g');
    ylabel('M_z [Nm]');

    % Gravity Gradient
    figure('Name',strcat(mission_name, ' Gravity Gradient Torques'));
    % M_x
    subplot(1,3,1);
    hold on;

    plot(downsampled_time, downsample(sim_output.ext_torques.grav_Mx, df));

    title_text = ['Gravity M_x of ', mission_name, ' in s/c principle axes'];
    title(title_text);
    xlabel(time_label);
    ytickformat('%.3g');
    ylabel('M_x [Nm]');

    % M_y
    subplot(1,3,2);
    hold on;

    plot(downsampled_time, downsample(sim_output.ext_torques.grav_My, df));

    title_text = ['Gravity M_y of ', mission_name, ' in s/c principle axes'];
    title(title_text);
    xlabel(time_label);
    ytickformat('%.3g');
    ylabel('M_y [Nm]');

    % M_z
    subplot(1,3,3);
    hold on;

    plot(downsampled_time, downsample(sim_output.ext_torques.grav_Mz, df));

    title_text = ['Gravity M_z of ', mission_name, ' in s/c principle axes'];
    title(title_text);
    xlabel(time_label);
    ytickformat('%.3g');
    ylabel('M_z [Nm]');

    % Magnetic
    figure('Name',strcat(mission_name, ' Magnetic Field'));
    % M_x
    subplot(1,3,1);
    hold on;

    plot(downsampled_time, downsample(sim_output.ext_torques.mag_Mx, df));

    title_text = ['Magnetic Field M_x of ', mission_name, ' in s/c principle axes'];
    title(title_text);
    xlabel(time_label);
    ytickformat('%.3g');
    ylabel('M_x [Nm]');

    % M_y
    subplot(1,3,2);
    hold on;

    plot(downsampled_time, downsample(sim_output.ext_torques.mag_My, df));

    title_text = ['Magnetic Field M_y of ', mission_name, ' in s/c principle axes'];
    title(title_text);
    xlabel(time_label);
    ytickformat('%.3g');
    ylabel('M_y [Nm]');

    % M_z
    subplot(1,3,3);
    hold on;

    plot(downsampled_time, downsample(sim_output.ext_torques.mag_Mz, df));

    title_text = ['Magnetic Field M_z of ', mission_name, ' in s/c principle axes'];
    title(title_text);
    xlabel(time_label);
    ytickformat('%.3g');
    ylabel('M_z [Nm]');
    
end
