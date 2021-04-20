
function plot_coordtriads(sim_constants, sim_output, plot_format)
    mission_name = plot_format.mission_name;
    N = 50; % only 100 triads per orbit, so axes are easy to see
    % Find last time index of first orbit using orbital E switch from 2pi to 0
    tend = find( diff(sign(sim_output.OE.E - pi)) == -2 , 1) - 1;
    df = fix(tend/N); 
    figure('Name',strcat(mission_name, ' Triads in ECI')); hold on;
    %% Plot Princ, Body Coord Triads in position in ECI
    subplot(1,2,1); hold on;
   
    % Get ECI locations
    X = downsample(sim_output.positions.x_123(1:tend,:), df);
    Y = downsample(sim_output.positions.y_123(1:tend,:), df);
    Z = downsample(sim_output.positions.z_123(1:tend,:), df);
    l = 0.9*norm([diff(X(1:2)) diff(Y(1:2)) diff(Z(1:2))]); % Length of coord axes in km
    % Make coordinate matrices (where columns are directions of vectors)
    R = downsample(sim_output.attitude.princ2inert(1:tend,:,:), df);
    W = zeros(size(R));
    for i = 1:size(R,1)
       W(i,:,:) = (sim_constants.rotm.')*squeeze(R(i,:,:)) ; % inertial->princ->body
    end

    % Plot inertial, principal, and body coordinates at each point
    plot_triad(X,Y,Z,R,l,[0 1 0],'Principal');
    plot_triad(X,Y,Z,W,l,[0 0 1],'Body');

    % Plot Earth for reference, but transparent
    R_Earth = sim_constants.R_Earth;
    [xE, yE, zE] = ellipsoid(0, 0, 0, R_Earth, R_Earth, R_Earth, 50);
    surface(xE, yE, zE, 'FaceColor', 'blue', 'EdgeColor', 'none','FaceAlpha',0.1,'DisplayName','Earth'); 

    title_text = ['Body & Principal Triads of ', mission_name, ...
        ' in ECI reference frame'];
    title(title_text);
    legend('Location','eastoutside');
    xlabel('X [km]');
    ylabel('Y [km]');
    zlabel('Z [km]');
    axis equal;
    pbaspect([1 1 1]);
    view(3);
    grid on;
    
    %% Plot RTN Coord Triad in position in ECI
    subplot(1,2,2); hold on;
    RTN = downsample(sim_output.positions.RTN2ECI(1:tend,:,:), df);
    plot_triad(X,Y,Z,RTN,l,[1 165/255 0],'RTN');
    % Plot Earth for reference, but transparent
    R_Earth = sim_constants.R_Earth;
    [xE, yE, zE] = ellipsoid(0, 0, 0, R_Earth, R_Earth, R_Earth, 50);
    surface(xE, yE, zE, 'FaceColor', 'blue', 'EdgeColor', 'none','FaceAlpha',0.1,'DisplayName','Earth'); 

    title_text = ['RTN Coordinate Triad of ', mission_name, ...
        ' in ECI'];
    title(title_text);
    legend('Location','eastoutside');
    xlabel('X [km]');
    ylabel('Y [km]');
    zlabel('Z [km]');
    axis equal;
    pbaspect([1 1 1]);
    view(3);
    grid on;
    
    %% Plot trace of coordinate tips for principal, body, and RTN axes in inertial frame
    figure('Name',strcat(mission_name, ' Triads in Time')); hold on;
    tend = fix(3*tend/100); % only plot a few percent of the orbital period to keep plot uncluttered
    sim_time = sim_output.time(1:tend)*plot_format.seconds_to_increment; 
    Nax = 2; % number of actual axes to plot (unit length arrows along XYZ)
    ax_reduction = fix(length(sim_time)/(Nax-1) - 1); % lower sampling for axes plotting
    time_df3 = downsample(sim_time,ax_reduction);
    axes_name = 'Principal Axes';
    R = sim_output.attitude.princ2inert(1:tend,:,:);
    h = [];
    g = [0;0;0];
    for i = 1:3
        h = [h; subplot(1,3,i)]; hold on;
        if i == 2
            for k = 1:length(R)
                R(k,:,:) = (sim_constants.rotm.')*squeeze(R(k,:,:)) ; % inertial->princ->body
            end
            axes_name = 'Body Axes';
        elseif i == 3
            R = sim_output.positions.RTN2ECI(1:tend,:,:);
            axes_name = 'RTN Axes';          
        end
        Xx = squeeze(R(:,1,1)); Yx = squeeze(R(:,1,2)); Zx = squeeze(R(:,1,3));
        Xy = squeeze(R(:,2,1)); Yy = squeeze(R(:,2,2)); Zy = squeeze(R(:,2,3));
        Xz = squeeze(R(:,3,1)); Yz = squeeze(R(:,3,2)); Zz = squeeze(R(:,3,3));
        timec = 1+round(255*sim_time/sim_time(end)); % convert to 255 mapping
        colormap spring;
        cmapx = colormap; 
        g(1) = patch([Xx; nan],[Xy; nan],[Xz; nan],reshape([cmapx(timec,:); [0 0 0]], [length(timec)+1,1,3]),'FaceColor','none','EdgeColor','interp','DisplayName','X Trace');
        colormap summer;
        cmapy = colormap; 
        g(2) = patch([Yx; nan],[Yy; nan],[Yz; nan],reshape([cmapy(timec,:); [0 0 0]], [length(timec)+1,1,3]),'FaceColor','none','EdgeColor','interp','DisplayName','Y Trace');
        colormap autumn;
        cmapz = colormap;
        g(3) = patch([Zx; nan],[Zy; nan],[Zz; nan],reshape([cmapz(timec,:); [0 0 0]], [length(timec)+1,1,3]),'FaceColor','none','EdgeColor','interp','DisplayName','Z Trace');
        Rdown = downsample(R, ax_reduction);
        timec = 1+round(255*time_df3/time_df3(end)); % convert to 255 mapping
        cmap = zeros(3, length(timec), 3);
        cmap(1,:,:) = cmapx(timec,:); cmap(2,:,:) = cmapy(timec,:); cmap(3,:,:) = cmapz(timec,:);
        plot_triad(zeros(length(Rdown),1),zeros(length(Rdown),1),zeros(length(Rdown),1),Rdown,1,cmap);
        xlabel('1');
        ylabel('2');
        zlabel('3');
        axis equal;
        pbaspect([1 1 1]);
        view(3);
        grid on;
        legend(g);
        title_text = [axes_name ' of ', mission_name, ...
            ' w.r.t. ECI'];
        title(title_text);
    end
%       cbh = colorbar(h(1)); 
%       h(1).Position(3:4) = h(2).Position(3:4);
%       set(get(cbh,'label'),'string','Time [sec]');
%      % Reposition to figure's left edge, centered vertically
%       cbh.Position(1) = .95-cbh.Position(3);
%       cbh.Position(2) = 0.5-cbh.Position(4)/2;

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
            c(c > 1) = 1;
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
            for j = 1:size(c,2)
                quiver3(X(j),Y(j),Z(j), xdir(j,1), xdir(j,2), xdir(j,3),'Color', squeeze(c(1,j,:)), ...
                            'LineWidth',wid,'AutoScale',aut,'DisplayName','');
                quiver3(X(j),Y(j),Z(j), ydir(j,1), ydir(j,2), ydir(j,3),'Color', squeeze(c(2,j,:)), ...
                           'LineWidth',wid,'AutoScale',aut,'DisplayName','');
                quiver3(X(j),Y(j),Z(j), zdir(j,1), zdir(j,2), zdir(j,3),'Color', squeeze(c(3,j,:)), ...
                           'LineWidth',wid,'AutoScale',aut,'DisplayName','');
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