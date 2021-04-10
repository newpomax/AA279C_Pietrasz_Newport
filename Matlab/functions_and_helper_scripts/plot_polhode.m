function plot_polhode(w0, I,wout,plot_name)
    N_faces = 500;
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
    subplot(3,3,[1 2 4 5 7 8]); hold on;
    [X1,X2,X3] = ellipsoid(0,0,0,energy_semiaxes(1),energy_semiaxes(2),energy_semiaxes(3),N_faces);
    h1 = surf(X1,X2,X3);
    set(h1,'FaceColor',[0 1 0],'FaceAlpha',0.5,'FaceLighting','gouraud','EdgeColor','none');
    % Plotting momentum ellipse
    [X1,X2,X3] = ellipsoid(0,0,0,mom_semiaxes(1),mom_semiaxes(2),mom_semiaxes(3),N_faces);
    h2 = surf(X1,X2,X3);
    set(h2,'FaceColor',[0 0 1],'FaceAlpha',0.5,'FaceLighting','gouraud','EdgeColor','none');
    % Plotting ohlode
    plot3(wout(:,1),wout(:,2),wout(:,3),'r','LineWidth',3);
    % Set legend
    legend('Energy Ellipse','Momentum Ellipse','Polhode');
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title('Polhode Plot');
    view(3); axis equal;
    
    %% Plot Pohlode along axes
    subplot(3,3,3); hold on;
    plot(wout(:,2),wout(:,3)); 
    xlabel('Y'); ylabel('Z'); axis equal;
    title('Pohlode Viewed Along X');
    subplot(3,3,6); hold on;
    plot(wout(:,1),wout(:,3)); 
    xlabel('X'); ylabel('Z'); axis equal;
    title('Pohlode Viewed Along Y');
    subplot(3,3,9); hold on;
    plot(wout(:,1),wout(:,2)); 
    xlabel('X'); ylabel('Y'); axis equal;
    title('Pohlode Viewed Along z');
end