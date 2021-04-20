% rot_body2pring
% Calculates principal moments of inertia and rotation matrix from body
% axes to principle axes (with increasing moment of inertia,Ix <= Iy <= Iz
%
% Input:
%       I - the 3x3 inertia matrix
%
% Output:
%       I_princ - a 3x1 vector of the principle moment of inertia,
%                  [Ixx, Iyy, Izz]
%       R - The rotation matrix from principle axes to body axes
function [I_princ, R] = rot_body2princ(I)
    [R, I_princ] = eig(I);
    I_princ = diag(I_princ);
    R_addl = eye(length(I_princ));
    if I_princ(1) > I_princ(2) % order the eigenvalues
        holdeig = I_princ(1);
        I_princ(1) = I_princ(2);
        I_princ(2) = holdeig;
        % we have flipped our x and y axes, so we need to add more rotation
        R_addl = rotz(-90)*rotx(180)*R_addl;
    end
    if I_princ(1) > I_princ(3)
        holdeig = I_princ(1);
        I_princ(1) = I_princ(3);
        I_princ(3) = holdeig;
        % we have flipped our x and z axes, so we need to add more rotation
        R_addl = roty(-90)*rotx(180)*R_addl;
    end
    if I_princ(2) > I_princ(3)
        holdeig = I_princ(2);
        I_princ(2) = I_princ(3);
        I_princ(3) = holdeig;
        % we have flipped our y and z axes, so we need to add more rotation
        R_addl = rotx(-90)*roty(180)*R_addl;
    end
    R = R/R_addl;
    R(:,2) = cross(R(3,:),R(1,:)); % ensure that trio is right-handed
end
    