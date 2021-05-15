%% EKF Testing Script
close all
clear all

% Add path to subfolders
addpath(fullfile('..', 'functions_and_helper_scripts'));

% Load simulation data
load('EKF_testing.mat');
w_meas = squeeze(w_meas.Data).';
B_meas = B_meas.Data;
B_meas = B_meas./(sqrt(sum(B_meas.^2,2))*ones(1,3)); % convert to just direction
S_meas = S_meas.Data;
B_ECI = B_ECI.Data;
B_ECI = B_ECI./(sqrt(sum(B_ECI.^2,2))*ones(1,3)); % convert to just direction
S_ECI = S_ECI.Data;
M = M_perturbations.Data; % torque input (should be p small generally)

% Add new sim_constants items
sim_constants.mu0 = [0;0;0;1;0.01;0.01;0.01;0;0;0;0;0]; % initial estimate
sim_constants.cov0 = eye(12); % initial covariance
sim_constants.Q = 0.1*sim_constants.time_step*eye(12); % process noise
sim_constants.R = 0.1*eye(9); % measurement noise (will need updated for new sensor measurements)
sim_constants.R(1:3,1:3) = diag(sim_constants.gyro_error);
sim_constants.R(4:6,4:6) = diag(sim_constants.mag_error);
% sim_constants.R = zeros(9);

% Run EKF simulator
sim('EKF_Testing');

%% Get estimated attitude in Euler312 angles
A_EKF = zeros(length(q.Data),3,3);
e_EKF = zeros(length(q.Data),3);
for i = 1:length(e_EKF)
    myA = q2A(q.Data(i,:).');
   A_EKF(i,:,:) = myA;
   e_EKF(i,:) = A2e(myA);
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

