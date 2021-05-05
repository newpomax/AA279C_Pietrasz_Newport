R_Earth = 6371;

function B_mag = magnetic_field(geod, r_ECI, R_Earth, gmst)
    % Calculates magnetic field, in ECI, at r, given in geod (deg, deg, km),
    % using a simple dipole model. % TODO: FIX ME :(
    
    lat = geod(1);
    colat = 90 - lat;
    long = geod(2);
    r = r(3);
    
    % Dipole properties
    a3H0 = 7.943*10^15; % Wb*m
    dipole_colat = 168.6; % deg
    dipole_long = 109.3; % deg
    dipole_a = 0; % TODO: calculate
    mag_Earth = [sind(dipole_colat)*cosd(dipole_a);
                 sind(dipole_colat)*sind(dipole_a);
                 cosd(dipole_colat)];
                  
    B = (a3H0/r^3)*(3*dot(mag_Earth,r_ECI)*r_ECI-mag_Earth);
    
    % Convert from geod to ECI
    RAAN = long + rad2deg(gmst);
    B_mag(1) = ( B_R*cosd(lat) + B_colat*cosd(lat) )*cosd(RAAN) ...
        - B_long*sind(RAAN);
    B_mag(2) = ( B_R*cosd(lat) + B_colat*cosd(lat) )*sind(RAAN) ...
        + B_long*sind(RAAN);
    B_mag(3) = B_R*sind(lat) - B_colat*cosd(lat);
    
end