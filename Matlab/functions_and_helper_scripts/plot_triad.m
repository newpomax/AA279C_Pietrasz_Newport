%% plot_triad
% plots a corrdinate trio defined by the rotation matrix R centered at 
% (x,y,z) in inertial space. l sets the length of the plotted coordinates
% (= 1 by default), c sets the base color for the first axis (the remaining
% two are slightly shifted to give distinction between axes).

function plot_triad(X, Y, Z, R, l, c, name)
    if nargin < 5
       l = 1; 
    end
    if nargin < 6
       c = [1 0 0]; 
    end
    if nargin < 7
        name = 'Coord Triad';
    end
    c = reshape(c,1,3);
    [~,maxi] = max(c);
    c = [c; abs(rot(maxi,15)*(c.')).';abs(rot(maxi,30)*(c.')).'];
    wid = 1;
    aut = 'off';
    % Make x
    xdir = l*squeeze(R(:,1,:))./sqrt(sum(squeeze(R(:,1,:)).^2,2)) ;
    % Make y
    ydir = l*squeeze(R(:,2,:))./sqrt(sum(squeeze(R(:,2,:)).^2,2)) ;
    % Make z
    zdir = l*squeeze(R(:,3,:))./sqrt(sum(squeeze(R(:,3,:)).^2,2)) ;
    quiver3(X,Y,Z, xdir(:,1), xdir(:,2), xdir(:,3),'Color', c(1,:), ...
                   'LineWidth',wid,'DisplayName',strcat(name, '_X'),'AutoScale',aut);
    quiver3(X,Y,Z, ydir(:,1), ydir(:,2), ydir(:,3),'Color', c(2,:), ...
                   'LineWidth',wid,'DisplayName',strcat(name, '_Y'),'AutoScale',aut);
    quiver3(X,Y,Z, zdir(:,1), zdir(:,2), zdir(:,3),'Color', c(3,:), ...
                   'LineWidth',wid,'DisplayName',strcat(name, '_Z'),'AutoScale',aut);
    function R = rot(ax, angle)
        switch ax(1)
            case 1
                R = roty(angle)*rotz(angle);
            case 2
                R = rotz(angle)*rotx(angle);
            case 3
                R = rotx(angle)*roty(angle);
        end
        R = abs(R);
    end
end