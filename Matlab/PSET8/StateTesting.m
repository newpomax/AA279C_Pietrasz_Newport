%% State testing
N = 1000;
S_inplane = ones(N,1)*[0 1 0];
angs = linspace(0,4*pi,N).';
v_ECI = [sin(angs) , cos(angs), zeros(N,1)];
w = zeros(N,3);
w(end-N/10:end,:) = ones(N/10+1,1)*[5 0 12];
w(end-N/20,:) = [5 0 5];

time_step = 1;
SLEW_LOWER = 0.05;
SLEW_UPPER = sind(30);
TUMBLE_LIMIT = 10; % rad/s
DEBOUNCE_TIME = 10*time_step; % state machine debounce time