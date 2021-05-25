%% Attitude Control Module Testing
clearvars; clc; close all;

N = 1000;
time_step = 1;
sim_constants.I_r = 1;
sim_constants.r_rotor = [1;0;0];
wr0 = 0;
dw_rdt = ones(N,1);
w_sc = ones(N,1)*deg2rad([0 0 3]);
q_sc = [ones(N/2,1)*[0 0 0 1] ; ones(N/2,1)*[0 sin(pi/4) 0 cos(pi/4)] ];
m_mag = ones(N,1)*[1 1 1]/sqrt(3);

sim('test_attcontrol');
M_mag = squeeze(M_mag.Data).';
M_r = squeeze(M_r.Data).';

figure('Name','M_mag'); hold on;
title('Magnetorquer unit testing');
tspan = 0:N-1;
labs = {'x','y','z'};
for i = 1:3
    subplot(1,3,i); hold on; title(['M_' labs{i}]);
    plot(tspan,M_mag(:,i));
    ylabel('Torque, Nm'); xlabel('Time');
end

figure('Name','M_r'); hold on;
title('Momentum wheel unit testing');
tspan = 0:N-1;
labs = {'x','y','z'};
for i = 1:3
    subplot(1,3,i); hold on; title(['M_' labs{i}]);
    plot(tspan,M_r(:,i));
    ylabel('Torque, Nm'); xlabel('Time');
end