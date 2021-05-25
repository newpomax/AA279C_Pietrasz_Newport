%% State testing
clearvars; close all; clc;
addpath(fullfile('..', 'functions_and_helper_scripts'));
N = 1000;
Aerr = zeros(N,3,3);
for i = 1:N/2
    Aerr(i,:,:) = eye(3) + deg2rad([0 3 -4;-3 0 2;4 -2 0]);
end
for i = N/2+1:N
    Aerr(i,:,:) = eye(3) + deg2rad([0 6 -4;-6 0 2;4 -2 0]);
end
w = zeros(N,3);
w(end-N/10:end,:) = ones(N/10+1,1)*[5 0 12];
w(end-N/20,:) = [5 0 5];

time_step = 1;
LARGE_ANGLE_LIM = deg2rad(5);
TUMBLE_LIMIT = 10; % rad/s
DEBOUNCE_TIME = 10*time_step; % state machine debounce time