disp('--------------------------------------------------');
disp('1st if, not 2nd if');
azy1 = 0.001
ely1 = pi/2-0.001
azg1 = 2*pi-0.001
elg1 = pi/2-0.001
[azy1, ely1, azg1, elg1] = z(azy1, ely1, azg1, elg1)

display('1st if, 2nd if, subif');
azy1 = pi/2-0.001
ely1 = pi/2-0.001
azg1 = 3*pi/2+0.001
elg1 = pi/2+0.001
[azy1, ely1, azg1, elg1] = z(azy1, ely1, azg1, elg1)

display('1st if, 2nd if, subelse');
azy1 = pi/2-0.001
ely1 = pi/2-0.001
azg1 = 3*pi/2+0.001
elg1 = -pi/2+0.001
[azy1, ely1, azg1, elg1] = z(azy1, ely1, azg1, elg1)

display('1st elif, not 2nd if');
azy1 = 2*pi-0.001
ely1 = pi/2-0.001
azg1 = 0.001
elg1 = pi/2-0.001
[azy1, ely1, azg1, elg1] = z(azy1, ely1, azg1, elg1)

display('1st elif, 2nd if, subif');
azy1 = 3*pi/2+0.001 % > 3pi/2, azy-azg > pi/2
ely1 = pi/2-0.001 
azg1 = pi/2-0.001 % < pi/2, azy-azg > pi/2
elg1 = pi/2+0.001 % > 0
[azy1, ely1, azg1, elg1] = z(azy1, ely1, azg1, elg1)

display('1st elif, 2nd if, subelse');
azy1 = 3*pi/2+0.001 % > 3pi/2, azy-azg > pi/2
ely1 = pi/2-0.001 
azg1 = pi/2-0.001 % < pi/2, azy-azg > pi/2
elg1 = -pi/2+0.001 % <= 0
[azy1, ely1, azg1, elg1] = z(azy1, ely1, azg1, elg1)

display('not 1st if/elif, 2nd if, subif');
azy1 = 0.001
ely1 = pi/2-0.001
azg1 = 3*pi/2-0.001
elg1 = pi/2+0.001
[azy1, ely1, azg1, elg1] = z(azy1, ely1, azg1, elg1)

display('not 1st if/elif, 2nd if, subelse');
azy1 = 0.001
ely1 = pi/2-0.001
azg1 = 3*pi/2-0.001
elg1 = -pi/2+0.001
[azy1, ely1, azg1, elg1] = z(azy1, ely1, azg1, elg1)
% note: ely < 0 does not happen bc of sensor model


function [azy, ely, azg, elg] = z(azy, ely, azg, elg)
    % check if azimuth is near the turn-over point 2*pi->0
    if azy < pi/2 && azg > 3*pi/2
       disp('1st if');
       azg = azg - 2*pi;
    elseif azy > 3*pi/2 && azg < pi/2
       disp('1st elif');
       azy = azy - 2*pi;
    end
    % check that azimuths aren't near-opposite, if so, flip
    if abs(azy-azg) > pi/2
        disp('2nd if');
        azg = mod(azg+pi, 2*pi); % flip azg
        % note: ely < 0 does not happen bc of sensor model
        if elg > 0
           disp('subif');
           elg = pi - elg;
        else
           disp('subelse');
           elg = -pi - elg;
        end
    end
    
    display('flipped! (or not)');
end