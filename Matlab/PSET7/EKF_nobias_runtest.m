%% EKF Testing Script
close all
clear all

% Add path to subfolders
addpath(fullfile('..', 'functions_and_helper_scripts'));

% Load simulation data
load('EKF_testing_nobias.mat');
w_meas = w_meas.Data;
B_meas = B_meas.Data;
B_meas = B_meas;%./(sqrt(sum(B_meas.^2,2))*ones(1,3)); % convert to just direction
Sazel_meas = Sazel_meas.Data;
B_ECI = B_ECI.Data;
B_ECI = B_ECI;%./(sqrt(sum(B_ECI.^2,2))*ones(1,3)); % convert to just direction
S_ECI = S_ECI.Data;
M = M_perturbations.Data; % torque input (should be p small generally)

% Add new sim_constants items
sim_constants.sunsensor_rotm = eye(3); % sun sensor frame = princ axes
sim_constants.mu0 = [0;0;0;1;deg2rad([0.1; 0.1; 3])]; % initial estimate
sim_constants.cov0 = eye(7); % initial covariance
sim_constants.Q = 1E-2*sim_constants.time_step*eye(7); % process noise
sim_constants.R = eye(8); % measurement noise 
sim_constants.R(1:3,1:3) = diag(sim_constants.gyro_error);
sim_constants.R(4:6,4:6) = diag(sim_constants.mag_error);
sim_constants.R(7:8,7:8) = diag(sim_constants.ss_error);
sim_constants.R = 1E2*sim_constants.R;

% Run EKF simulator
sim('EKF_noBias_Testing');

%% Get estimated attitude in Euler312 angles
A_EKF = zeros(length(q_EKF.Data),3,3);
e_EKF = zeros(length(q_EKF.Data),3);
e_EKFesterr = zeros(size(e_EKF));
e_EKFtargerr = zeros(size(e_EKF));
q_cov = zeros(length(q_EKF.Data),4);
w_cov = zeros(size(e_EKF));
for i = 1:length(e_EKF)
    myq = q_EKF.Data(i,:);
    myq = myq/norm(myq);
    myA = q2A(myq.');
   A_EKF(i,:,:) = myA;
   e_EKFesterr(i,:) = A2e(squeeze(A.Data(:,:,i))*(myA.'));
   e_EKFtargerr(i,:) = A2e(squeeze(A_target.Data(:,:,i))*(myA.'));
   
   e_EKF(i,:) = A2e(myA);
   q_cov(i,:) = diag(squeeze(cov_q.Data(:,:,i)));
   w_cov(i,:) = diag(squeeze(cov_w.Data(:,:,i)));
end

%% Plot stuff
time_label = ['Time [', num2str(plot_format.time_increments), ']'];
df = plot_format.downsample_factor;
time_df = downsample(sim_output.time,df);
figure('Name','Attitude Estimate'); hold on;
labs = {'\phi', '\theta', '\psi'};

for i = 1:3
   subplot(1,3,i); hold on;
   plot(time_df,wrapTo180(rad2deg(downsample(e_estq.Data(:,i),df))), ':', 'Linewidth',1,'DisplayName','q method');
   plot(time_df,wrapTo180(rad2deg(downsample(e.Data(:,i),df))),'LineWidth',1,'DisplayName','Attitude Angle');
   plot(time_df,wrapTo180(rad2deg(downsample(e_EKF(:,i),df))),':','LineWidth',1.5,'DisplayName','EKF');
   legend; ylabel([labs{i} ', deg']); xlabel(time_label);
   title_text = [labs{i} ' of ', plot_format.mission_name];
   title(title_text); legend('location','best');
end

figure('Name','Attitude Est Err'); hold on;
labs = {'\phi', '\theta', '\psi'};

for i = 1:3
   subplot(1,3,i); hold on;
   plot(time_df,wrapTo180(rad2deg(downsample(e_esterr.Data(:,i),df))),'LineWidth',1,'DisplayName','Est. Err, Q-method');
   plot(time_df,wrapTo180(rad2deg(downsample(e_EKFesterr(:,i),df))),':','LineWidth',1.5,'DisplayName','ESt. Err, EKF');
   legend; ylabel([labs{i} ', deg']); xlabel(time_label);
   title_text = [labs{i} ' of ', plot_format.mission_name];
   title(title_text); legend('location','best');
end

figure('Name','Ang Vel Estimate'); hold on;
labs = {'\omega_x', '\omega_y', '\omega_z'};
for i = 1:3
   subplot(1,3,i); hold on;
   plot(time_df,rad2deg(downsample(w_meas(:,i),df)), ':', 'Linewidth',1,'DisplayName','Measured \omega');
   plot(time_df,rad2deg(downsample(w.Data(:,i),df)),'LineWidth',1,'DisplayName','True \omega');
   plot(time_df,rad2deg(downsample(w_EKF.Data(:,i),df)),':','LineWidth',1.5,'DisplayName','EKF');
   legend; ylabel([labs{i} ', deg/s']); xlabel(time_label);
   title_text = [labs{i} ' of ', plot_format.mission_name];
   title(title_text); legend('location','best');
end

figure('Name','Ang Vel Resids'); hold on;
labs = {'\omega_1', '\omega_2', '\omega_3'};
for i = 1:3
   subplot(1,3,i); hold on;
   plot(time_df,downsample(z_pre.Data(:,i),df), ':', 'Linewidth',1,'DisplayName','Pre-fit Difference');
   plot(time_df,downsample(z_post.Data(:,i),df), ':', 'Linewidth',1,'DisplayName','Post-fit Difference');
   legend; ylabel([labs{i}]); xlabel(time_label);
   title_text = [labs{i} ' residuals of ', plot_format.mission_name];
   title(title_text); legend('location','best');
end

figure('Name','Mag Resids'); hold on;
labs = {'B_1', 'B_2', 'B_3'};
for i = 1:3
   subplot(1,3,i); hold on;
   plot(time_df,downsample(z_pre.Data(:,i+3),df), ':', 'Linewidth',1,'DisplayName','Pre-fit Difference');
   plot(time_df,downsample(z_post.Data(:,i+3),df), ':', 'Linewidth',1,'DisplayName','Post-fit Difference');
   legend; ylabel([labs{i} ', Tesla']); xlabel(time_label);
   title_text = [labs{i} ' residuals of ', plot_format.mission_name];
   title(title_text); legend('location','best');
end

figure('Name','SS Resids'); hold on;
labs = {'azimuth','elevation'};
for i = 1:2
   subplot(1,2,i); hold on;
   plot(time_df,rad2deg(downsample(z_pre.Data(:,i+6),df)), ':', 'Linewidth',1,'DisplayName','Pre-fit Difference');
   plot(time_df,rad2deg(downsample(z_post.Data(:,i+6),df)), ':', 'Linewidth',1,'DisplayName','Post-fit Difference');
   legend; ylabel([labs{i} ', deg']); xlabel(time_label);
   title_text = [labs{i} ' residuals of ', plot_format.mission_name];
   title(title_text); legend('location','best');
end