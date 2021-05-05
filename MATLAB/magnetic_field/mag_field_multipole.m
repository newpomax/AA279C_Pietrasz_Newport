R_Earth = 6371;
igrf = readmatrix('igrf13coeffs.txt');
igrf = igrf(2:25,28); % IGRF 2020 model coefficients
% First N=4 coefficients ordered by n, then m, then g/h

% B_mag = magnetic_field([0, 0, R_Earth], R_Earth, igrf, 0)
% norm(B_mag)
B_mag = magnetic_field([90, 0, R_Earth], R_Earth, igrf, 0)
norm(B_mag)
B_mag = magnetic_field([0, 90, R_Earth], R_Earth, igrf, 0)
norm(B_mag)
% B_mag = magnetic_field([0, 0, R_Earth+100], R_Earth, igrf, 0)
% norm(B_mag)
% B_mag = magnetic_field([0, 0, R_Earth+100], R_Earth, igrf, 180)
% norm(B_mag)

function B_mag = magnetic_field(r, R_Earth, igrf, gmst)
    % Calculates magnetic field, in ECI, at r, given in geod (deg, deg, km),
    % using a 4th-order pole model. 
    % https://www.ngdc.noaa.gov/IAGA/vmod/igrf.html
    
    k_max = 4;
    
    lat = r(1);
    colat = 90 - lat;
    long = r(2);
    r = r(3);
    
    % Calculate potential for a small region around r in order to take the
    % numerical gradient
    d = 0.01;
    colat_vec = [colat-d:d:colat+d];
    long_vec = [long-d:d:long+d];
    r_vec = [r-100*d:100*d:r+100*d];
    
    % Magnetic potential as mesh of spherical harmonics
    V = zeros(3,3,3);

    for n = 1:k_max
        V_n = zeros(3,3,3);
        
        % Extract relevant coefficients
        gh_coeff = igrf(n^2:(n+1)^2-1);
        g = gh_coeff(1:2:end);
        h = gh_coeff(2:2:end);
        
        % Calculate Legendre function
        V_m = zeros(3,3);
        P_mn = legendre(n, cosd(colat_vec));
        
        % Multiply coefficients
        for m=0:n
            if m == 0
                % Transpose .' required to correctly compute outer product
                V_m = V_m + g(m+1)*cosd(long_vec).' * P_mn((m+1),:);
            else
                V_m = V_m + ( g(m+1)*cosd((m+1)*long_vec) + ...
                    h(m)*sind((m+1)*long_vec) ).' * P_mn(m+1,:);
            end
        end
        rn_normalized = (R_Earth./r_vec).^(n+1);
        % Outer product of rn_normalized and V_m
        for i = 1:3
            V_n(i,:,:) = rn_normalized(i)*V_m;
        end
        
        V = V + V_n;
    end
    V = V*R_Earth; % in nT
    
    % Find field from potential. Use (2,2,2) as the center of our gradient,
    % Transform from partial derivatives of V wrt geod to grad V wrt geod. 
    [V_R, V_colat, V_long] = gradient(V);
    B_colat = -(1/r)*V_colat(2,2,2);
    B_long = -(1/r*sin(colat))*V_long(2,2,2);
    B_R = -V_R(2,2,2);
    
    % Convert from geod to ECI
    RAAN = long + rad2deg(gmst);
    B_mag(1) = ( B_R*cosd(lat) + B_colat*cosd(lat) )*cosd(RAAN) ...
        - B_long*sind(RAAN);
    B_mag(2) = ( B_R*cosd(lat) + B_colat*cosd(lat) )*sind(RAAN) ...
        + B_long*sind(RAAN);
    B_mag(3) = B_R*sind(lat) - B_colat*cosd(lat);
    
end