function PlotWing(dx,c,b,Ny,Nx,y_panel,x_panel,xsc0)

x_LE_sc = dx;        % LE = dx(y)
x_TE_sc = c + dx;    % TE = c(y) + dx(y)

figure; hold on; axis equal; box on;

patch([b, fliplr(b)],[x_LE_sc, fliplr(x_TE_sc)],zeros(1,2*(Ny+1)), ...
      'FaceColor','none','EdgeColor','k');
% 1/4 chord lines
for i = 1:Nx
    eta_q = (i-3/4)/Nx;
    x_chord_sc = dx + eta_q*c;
    plot(b, x_chord_sc, 'b-');
end
% Colocation points
plot(y_panel(:,3),x_panel(:,3),'.r');
% Shear center line  (constant)
plot(b, xsc0*ones(size(b)), 'k-', 'LineWidth', 2);

% Spanwise lines of elements, from root to tip
for i = 0:Nx
    eta = i/Nx;
    x_line = dx+ eta*c;
    plot(b,x_line,'--','Color',[0.4 0.4 0.4],'LineWidth', 0.5);
end

% Span divisions, vertical lines
for j = 1:(Ny+1)
    x_seg = linspace(x_LE_sc(j), x_TE_sc(j), 50);  % de LE a TE
    y_seg = b(j)*ones(size(x_seg));             % span fijo
    plot(y_seg,x_seg,'--','Color',[0.4 0.4 0.4],'LineWidth', 0.5);
end

ylim([-1 1.5]);
xlabel('Spanwise y');
ylabel('x');
title('Mesh with shear center aligned');
end