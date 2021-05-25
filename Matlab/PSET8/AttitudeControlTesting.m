%% Attitude Control Module Testing
clearvars; clc; close all;

sim_constants.dw_rdt_err = (deg2rad(0.1)/2).^2; %rad/s^2
sim_constants.dw_rdt_bias = deg2rad(0.1); %rad/s^2
sim_constants.m_magtor_err = ([.01 .01 .01]/2).^2; % A*m^2
sim_constants.m_magtor_bias = [-.01 .005 -.02]; % A*m^2

N = 1000;
time_step = 1;
sim_constants.I_r = 0.01;
sim_constants.r_rotor = [1;0;0];
wr0 = 0;
dw_rdt = ones(N,1);
w_sc = ones(N,1)*deg2rad([0 0 1]);
q_sc = [ones(N/2,1)*[0 0 0 1] ; ones(N/2,1)*[0 sin(pi/4) 0 cos(pi/4)] ];
m_mag = ones(N,1)*[1 1 1]/sqrt(3);

sim('test_attcontrol');
M_mag = squeeze(M_mag.Data).';
M_mag_noerr = squeeze(M_mag_noerr.Data).';
M_r = squeeze(M_r.Data).';
M_r_noerr = squeeze(M_r_noerr.Data).';

figure('Name','M_mag'); hold on;
title('Magnetorquer unit testing');
tspan = 0:N-1;
labs = {'x','y','z'};
for i = 1:3
    subplot(1,3,i); hold on; title(['M_' labs{i}]);
    plot(tspan,M_mag(:,i));
    plot(tspan,M_mag_noerr(:,i));
    ylabel('Torque, Nm'); xlabel('Time');
    legend('True Torque','Expected Torque');
end

figure('Name','M_r'); hold on;
title('Momentum wheel unit testing');
tspan = 0:N-1;
labs = {'x','y','z'};
for i = 1:3
    subplot(1,3,i); hold on; title(['M_' labs{i}]);
    plot(tspan,M_r(:,i));
    plot(tspan,M_r_noerr(:,i));
    ylabel('Torque, Nm'); xlabel('Time');
    legend('True Torque','Expected Torque');
end