%% target_attitude
% Calculates direction cosine matrix from desired principal axes to
% inertial axes (ECI). 
function A_target = target_attitude(r_ECI, v_ECI, S, princ2body)
    % Calculate in-plane sun direction
    h = cross(r_ECI, v_ECI); % angular momentum vector, normal to orbital plane
    h = h/norm(h); % unitize
    inplane_sun_dir = S - dot(h,S)*h; % get in-plane sun direction
    if all(abs(inplane_sun_dir) < 1E-6)% if sun direction is perpendicular to orbital plane 
        inplane_sun_dir = r_ECI/norm(r_ECI); % just assume sun is in radial direction
    else
        inplane_sun_dir = inplane_sun_dir/norm(inplane_sun_dir);
    end
    A_body2ECI = zeros(3);
    if dot(v_ECI, S) < 0 % if moving away from the sun
        A_body2ECI(:,2) = h; % align body y axis with orbit plane normal
        A_body2ECI(:,3) = -inplane_sun_dir; % align body z to be along projection of -S in orbit plane (SRP thrusting in-plane)
        A_body2ECI(:,1) = cross( A_body2ECI(:,2) ,  A_body2ECI(:,3) ); % x completes triad
    else
        A_body2ECI(:,2) = h; % align body y axis with orbit plane normal
        A_body2ECI(:,1) = inplane_sun_dir; % align body x to be along projection of S in orbit plane (no thrusting)
        A_body2ECI(:,3) = cross( A_body2ECI(:,1) ,  A_body2ECI(:,2) ); % z completes triad
    end
    A_target = A_body2ECI*(princ2body.'); % princ->body->ECI
end