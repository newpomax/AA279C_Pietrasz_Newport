function plot_polhode(w0, I,wout,plot_name)
    N_faces = 1000;
    if nargin < 4
        figure('Name','Polhode'); hold on;
    else
        figure('Name',plot_name); hold on;
    end
    L = sqrt(sum((I.*w0).^2)); % angular momentum, magnitude 
    twoT = sum(I.*(w0.*w0)); % angular rotational energy
    energy_semiaxes = sqrt(twoT./I);
    mom_semiaxes = L./(I);
    
    %% Plot 3D ellipses & polhode
    % Plotting energy ellipse
    subplot(3,5,[1 2 3 6 7 8 11 12 13]); hold on;
    [X1,X2,X3] = ellipsoid(0,0,0,energy_semiaxes(1),energy_semiaxes(2),energy_semiaxes(3),N_faces);
    h1 = surf(X1,X2,X3);
    set(h1,'FaceColor',[0 1 0],'FaceAlpha',0.5,'FaceLighting','gouraud','EdgeColor','none');
    % Plotting momentum ellipse
    [X1,X2,X3] = ellipsoid(0,0,0,mom_semiaxes(1),mom_semiaxes(2),mom_semiaxes(3),N_faces);
    h2 = surf(X1,X2,X3);
    set(h2,'FaceColor',[0 0 1],'FaceAlpha',0.5,'FaceLighting','gouraud','EdgeColor','none');
    % Plotting ohlode
    plot3(wout(:,1),wout(:,2),wout(:,3),'r-o','LineWidth',2,'MarkerSize',2);

    % Set legend
    legend('Energy Ellipse','Momentum Ellipse','Polhode');
    xlabel('\omega_x, rad/s'); ylabel('\omega_y, rad/s'); zlabel('\omega_z, rad/s');
    title(sprintf('Polhode Plot, Initial Ang Vel [%0.1f, %0.1f, %0.1f] deg/s', rad2deg(w0)));
    view(3); axis equal;
    
    %% Plot Pohlode along axes
    subplot(3,5,[4 5]); hold on;
    plot(wout(:,2),wout(:,3),'-o','MarkerSize',2); 
    xlabel('\omega_y, rad/s'); ylabel('\omega_z, rad/s'); axis equal;
    title('Pohlode Viewed Along X');
    subplot(3,5,[9 10]); hold on;
    plot(wout(:,1),wout(:,3),'-o','MarkerSize',2); 
    xlabel('\omega_x, rad/s'); ylabel('\omega_z, rad/s'); axis equal;
    title('Pohlode Viewed Along Y');
    subplot(3,5,[14 15]); hold on;
    plot(wout(:,1),wout(:,2),'-o','MarkerSize',2); 
    xlabel('\omega_x, rad/s'); ylabel('\omega_y, rad/s'); axis equal;
    title('Pohlode Viewed Along z');
end