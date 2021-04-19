
function plot_coordtriads(sim_constants, sim_output, plot_format)
    mission_name = plot_format.mission_name;
    time_label = ['Time [', num2str(plot_format.time_increments), ']'];
    N = 100; % only 100 triads per orbit, so axes are easy to see
    df = fix(length(sim_output.time)/N); 
    % Find last time index of first orbit using orbital E switch from 2pi to 0
    tend = find( diff(sign(sim_output.OE.E - pi)) == -2 , 1) - df;
    figure('Name',strcat(mission_name, ' Coordinate Triads')); hold on;
    %% Plot Coord Triads in position in ECI
    subplot(2,3,[1 2 4 5]); hold on;
   
    % Get ECI locations
    X = downsample(sim_output.positions.x_123(1:tend,:), df);
    Y = downsample(sim_output.positions.y_123(1:tend,:), df);
    Z = downsample(sim_output.positions.z_123(1:tend,:), df);
    l = 0.9*norm([diff(X(1:2)) diff(Y(1:2)) diff(Z(1:2))]); %500; % Length of coord axes in km
    % Make coordinate matrices (where columns are directions of vectors)
    R = downsample(sim_output.attitude.princ2inert(1:tend,:,:), df);
    W = zeros(size(R));
    Q = zeros(size(R));
    for i = 1:size(R,1)
       W(i,:,:) = (sim_constants.rotm.')*squeeze(R(i,:,:)) ; % inertial->princ->body
       Q(i,:,:) = eye(3); % inertial -> inertial
    end

    % Plot inertial, principal, and body coordinates at each point
    plot_triad(X,Y,Z,Q,l,[1 0 0],'Inertial');
    plot_triad(X,Y,Z,R,l,[0 1 0],'Principal');
    plot_triad(X,Y,Z,W,l,[0 0 1],'Body');
    % Plot Earth for reference, but transparent
    R_Earth = sim_constants.R_Earth;
    [xE, yE, zE] = ellipsoid(0, 0, 0, R_Earth, R_Earth, R_Earth, 50);
    surface(xE, yE, zE, 'FaceColor', 'blue', 'EdgeColor', 'none','FaceAlpha',0.1,'DisplayName','Earth'); 

    title_text = ['Coordinate Triads of ', mission_name, ...
        ' in ECI reference frame'];
    title(title_text);
    legend('Location','westoutside');
    xlabel('X [km]');
    ylabel('Y [km]');
    zlabel('Z [km]');
    axis equal;
    pbaspect([1 1 1]);
    view(3);
    grid on;
    
    %% Plot trace of coordinate tips for principal, body axes in inertial frame
    tend = fix(3*tend/100); % only plot a few percent of the orbital period to keep plot uncluttered
    axes_name = 'Principal Axes';
    sim_time = sim_output.time(1:tend); 
    Nax = 3; % number of actual axes to plot (unit length arrows along XYZ)
    ax_reduction = round(length(sim_time)/Nax); % lower sampling for axes plotting
    time_df3 = downsample(sim_time, ax_reduction); 
    for i = 1:2
        subplot(2,3,3*i); hold on;
        R = sim_output.attitude.princ2inert(1:tend,:,:);
        if i == 2
            for k = 1:length(R)
                R(k,:,:) = (sim_constants.rotm.')*squeeze(R(k,:,:)) ; % inertial->princ->body
            end
            axes_name = 'Body Axes';
        end
        Xx = squeeze(R(:,1,1)); Yx = squeeze(R(:,1,2)); Zx = squeeze(R(:,1,3));
        Xy = squeeze(R(:,2,1)); Yy = squeeze(R(:,2,2)); Zy = squeeze(R(:,2,3));
        Xz = squeeze(R(:,3,1)); Yz = squeeze(R(:,3,2)); Zz = squeeze(R(:,3,3));
        patch([Xx; nan],[Xy; nan],[Xz; nan],[sim_time; nan],'FaceColor','none','EdgeColor','interp','DisplayName','X Trace');
        patch([Yx; nan],[Yy; nan],[Yz; nan],[sim_time; nan],'FaceColor','none','EdgeColor','interp','DisplayName','Y Trace');
        patch([Zx; nan],[Zy; nan],[Zz; nan],[sim_time; nan],'FaceColor','none','EdgeColor','interp','DisplayName','Z Trace');
        h = colorbar('westoutside');

        h.Position = h.Position - [0.08 0 0 0];
        set(get(h,'label'),'string',time_label);
        R = downsample(R, ax_reduction);
        cmap = colormap;
        plot_triad(zeros(length(R),1),zeros(length(R),1),zeros(length(R),1),R,1,cmap(1+round(255*time_df3/time_df3(end)),:));
        xlabel('1');
        ylabel('2');
        zlabel('3');
        axis equal;
        pbaspect([1 1 1]);
        view(3);
        grid on;
        title_text = [axes_name ' of ', mission_name, ...
            ' in inertial reference frame'];
        title(title_text);
    end
    %% plot_triad
    % plots a coordinate trio defined by the rotation matrix R centered at 
    % (x,y,z) in inertial space. l sets the length of the plotted coordinates
    % (= 1 by default), c sets the base color for the first axis (the remaining
    % two are slightly shifted to give distinction between axes).
    function plot_triad(X, Y, Z, R, l, c, name)
        rotc = false;
        if size(c,1) == 1 % if only provided one color, rotate axis colors to be slightly different
            rotc = true;
            [~,maxi] = max(c);
            c = [c; abs(rot(maxi,15)*(c.')).'; abs(rot(maxi,30)*(c.')).'];
        end
        wid = 1;
        aut = 'off';
        % Make x
        xdir = l*squeeze(R(:,:,1)) ;
        % Make y
        ydir = l*squeeze(R(:,:,2)) ;
        % Make z
        zdir = l*squeeze(R(:,:,3)) ;
        if rotc
            quiver3(X,Y,Z, xdir(:,1), xdir(:,2), xdir(:,3),'Color', squeeze(c(1,:)), ...
                       'LineWidth',wid,'DisplayName',strcat(name, '_X'),'AutoScale',aut);
            quiver3(X,Y,Z, ydir(:,1), ydir(:,2), ydir(:,3),'Color', squeeze(c(2,:)), ...
                       'LineWidth',wid,'DisplayName',strcat(name, '_Y'),'AutoScale',aut);
            quiver3(X,Y,Z, zdir(:,1), zdir(:,2), zdir(:,3),'Color', squeeze(c(3,:)), ...
                       'LineWidth',wid,'DisplayName',strcat(name, '_Z'),'AutoScale',aut);
        else
            for j = 1:length(c)
                quiver3(X(j),Y(j),Z(j), xdir(j,1), xdir(j,2), xdir(j,3),'Color', squeeze(c(j,:)), ...
                            'LineWidth',wid,'AutoScale',aut);
                quiver3(X(j),Y(j),Z(j), ydir(j,1), ydir(j,2), ydir(j,3),'Color', squeeze(c(j,:)), ...
                           'LineWidth',wid,'AutoScale',aut);
                quiver3(X(j),Y(j),Z(j), zdir(j,1), zdir(j,2), zdir(j,3),'Color', squeeze(c(j,:)), ...
                           'LineWidth',wid,'AutoScale',aut);
            end
        end
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
end