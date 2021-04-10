%% Testing EulerEq equations
close all; clearvars; clc;

%% Constants
% Add path to subfolders
addpath(fullfile('..', 'functions_and_helper_scripts'));

% Load constants
constants;

%% Simulation
output = sim('EulerEqs',1000);
plot_polhode(sim_constants.angvel0,sim_constants.I_princ,output.w);

%% Try a different angular vel
old_w0 = sim_constants.angvel0;
sim_constants.angvel0 = deg2rad([-2; .2; 8]); %rad/s, initial angular velocity vector in princ axes
output = sim('EulerEqs',100);
plot_polhode(sim_constants.angvel0,sim_constants.I_princ,output.w, 'Other Ang Vel');

sim_constants.angvel0 = old_w0; % return constant to its old value