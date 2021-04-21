%% OE2RTN -
% Takes in orbital elements and spits out the rotation matrix from RTN frame
% to ECI. 
function R = OE2RTN(a,e,i,RAAN,w,M,mu_g)
    R = zeros(3,3);
    if e == 0 % protect against circular orbit (undefined w)
       w = 0; 
    end
    if i == 0 % protect against equatorial orbit (undefined i)
       RAAN = 0; 
    end
    %Calculate additional useful orbit parameters
    E = mean2eccentric(M,e,1E-8);
    p = a*(1-e^2);
    h = sqrt(mu_g * p);
    
    x_PQW = a*(cos(E) - e);
    y_PQW = a*sqrt(1-e^2)*sin(E);
    z_PQW = 0;
    
    %Calculate position and velocity vectors
    rotm = rotz(-RAAN)*rotx(-i)*rotx(-w);
        % ... Rotation matrix from PQW to IJK 
    hECI = (rotm*([0;0;h])); % rotate ang momentum from PQW to IJK
    rECI = rotm*[x_PQW;y_PQW;z_PQW]; % rotate ECI from PQW to IJK
    R(:,1) = rECI/norm(rECI); % normalize radial vector
    R(:,3) = hECI/norm(hECI); % parallel to angular momentum
    R(:,2) = cross(R(:,3),R(:,1)); % T completes the right-handed triad
    
    function E = mean2eccentric(M, e, tolerance)
        % Given the mean anomaly M (in rad), eccentricity e, and tolerance, 
        % calculates the eccentric anomaly E (in rad) using the Newton-Raphson
        % method.

        % Initialize N-R
        E = M;
        d = calculate_new_delta(E, M, e);

        while abs(d) >= tolerance
            E = E + d;
            d = calculate_new_delta(E, M, e);
        end
    end

    function d = calculate_new_delta(E, M, e)
        % Calculates a new delta d to be used for calculating the eccentric 
        % anomaly E (in rad) using the Newton-Raphson method.

        d = -(E - e*sin(E) - M)/(1 - e*cos(E));
    end
end