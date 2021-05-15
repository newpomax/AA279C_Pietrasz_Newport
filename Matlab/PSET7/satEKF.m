%% satEKF
% Extended Kalman filter for satellite attitude state estimation.
% mu_in : the current estimated state mean, in the form
%           (q, w, bB, bS)^T
% cov_in : the current covariance matrix
% torque : the input torque to the system at the current time step
% meas : the measurement of the current state
% B_ECI : the modeled magnetic field vector in inertial coordinates at current time
% S_ECI: the modeled sun-pointing vector in intertial coordinates at current time 
% Q : the process noise matrix
% R : the measurement noise matrix
% I_princ : the principal moments of inertia corresponding to principal XYZ
% dt : the time step
function [mu_out, cov_out, z_pre, z_post] = satEKF(mu_in, cov_in, torque, meas, B_ECI, S_ECI, Q, R, I_princ, dt)
    q = mu_in(1:4); % quaternion
    w = mu_in(5:7); % angular velocity in princ
%     bB = mu_in(8:10); % component-wise bias of magnetic measurements 
    bS = mu_in(11:12); % azimuth and elevation bias of Sun sensor measurements
    
    At = A(q,w,I_princ);
    
    % Compute small-time quaternion change from current time to update time
    wnorm = norm(w);
    qnew = (eye(4)*cos(wnorm*dt*0.5) + Omega(w)*sin(wnorm*dt*0.5)/wnorm)*q;
    
    % Get approximate rotation matrix from ECI->princ at update time
    Anew = q2A(qnew); 
    %% Predict - use non-linear function to predict mean, linearized info for cov
    mu_p = f(mu_in, torque, I_princ);
    cov_p = At*cov_in*(At.') + Q;
    meas_p = g(mu_p,Anew*B_ECI,Anew*S_ECI); % predicted measurement
    z_pre = meas - meas_p; % pre-fit residual (diff between real and predicted)
    %% Update - use predicted states to compute Kalman gain and final update
    
    % Update state using sensitivity matrix
    Ct = C(q,bS,Anew,B_ECI,S_ECI); % Sensitivity matrix
    K = cov_p*(Ct.')/(Ct*cov_p*(Ct.') + R) ; % Kalman Gain
    mu_out = mu_p + K*z_pre;
    cov_out = cov_p - K*Ct*cov_p ;
    z_post = meas - g(mu_out,Anew*B_ECI,Anew*S_ECI); % post-fit residual (diff between real and estimated)
end
% Helper Functions
    function O = Omega(w)
        % Returns Omega(w), the rate relation such that dw/dt = 0.5*Omega*q
        O = [0 w(3) -w(2) w(1);
             -w(3) 0 w(1) w(2);
             w(2) -w(1) 0 w(3);
             -w(1) -w(2) -w(3) 0];
    end
    function L = Lambda(q)
        % Returns Lambda(q), the rate relation such that dw/dt = 0.5*Lambda*w
        L = [q(4) -q(3) q(2);
             q(3) q(4) -q(1);
             -q(2) q(1) q(4);
             -q(1) -q(2) -q(3) ];
    end
    function D = Delta(w,I)
        % Returns Delta(w), the partial derivative of the Euler equations
        % wrt angular velocity.
        D = [0 (I(2)-I(3))/I(1)*w(3) (I(2)-I(3))/I(1)*w(2);
             (I(3)-I(1))/I(2)*w(3) 0 (I(3)-I(1))/I(2)*w(1);
             (I(1)-I(2))/I(3)*w(2) (I(1)-I(2))/I(3)*w(1) 0 ];
    end
    function G = Gamma(q,V)
        % Returns Gamma(q,V), the partial derivative of A(q)*V wrt to the
        % quaternion, for some vector V in ECI.
        qvec = q(1:3);
        a = [-q(2);q(1);-q(4)];
        b = [-q(3);q(4);q(1)];
        c = [q(4); q(3); -q(2)];
        G = 2*[dot(qvec,V) dot(a,V) dot(b,V) dot(c,V);
               dot(-a,V) dot(qvec,V) dot(-c,V) dot(b,V);
               dot(-b,V) dot(c,V) dot(qvec,V) dot(-a,V) ];
    end
    function Mt = M(S)
        S12 = norm(S(1:2));
        Mt = [-S(2) S(1)*S(3)/S12;
               S(1) S(2)*S(3)/S12;
               0    -S12 ];
    end
    function Mt = Mtilde(b,S)
        S12 = norm(S(1:2));
        Mt = [b(2)*(S(2)^2)*S(3)/(S12^3) -b(1)-b(2)*S(1)*S(2)*S(3)/(S12^3) b(2)*S(1)/S12;
              b(1)-b(2)*S(1)*S(2)*S(3)/(S12^3) b(2)*(S(1)^2)*S(3)/(S12^3)  b(2)*S(2)/S12;
              -b(2)*S(1)/S12 -b(2)*S(2)/S12 0 ];
    end
    function xout = f(x,torque,I)
        % Dynamics function such that x_{t+1} = f(x_t, u_t)
        q = x(1:4);
        w = x(5:7);
        xout = [0.5*Omega(w)*q ;
               (torque - cross(w,I.*w))./I ;
               zeros(5,1) ];
    end
    function At = A(q,w,I)
        % Jacobian of dynamics
        At = [0.5*Omega(w) 0.5*Lambda(q) zeros(4,5);
              zeros(3,4) Delta(w,I) zeros(3,5);
              zeros(5,12) ] ;
    end
    function y = g(x,B_XYZ,S_XYZ)
        % Measurement function such that y_t = g(x_t,u_t)
        % A is the computed rotation matrix taken from the first four
        % elements of x
        y = zeros(9,1);
        y(1:3) = x(5:7); % angular velocity measurement
        y(4:6) = B_XYZ; % magnetic field strength in XYZ
        y(7:9) = S_XYZ; % sun direction in XYZ
    end
    function Ct = C(q,bS,Aq,B_ECI,S_ECI)
        S_xyz = Aq*S_ECI;
        Ct = [zeros(3,4) eye(3) zeros(3,5) ;
              Gamma(q,B_ECI) zeros(3) eye(3) zeros(3,2);
              (eye(3)+Mtilde(bS,S_xyz))*Gamma(q,S_ECI) zeros(3,6) M(S_xyz)];
    end