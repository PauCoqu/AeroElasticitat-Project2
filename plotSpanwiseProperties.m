function plotSpanwiseProperties(span)
y = span.y(:);

figure('Name','Cross-section properties along the span','Color','w');

tiledlayout(3,3,'Padding','compact','TileSpacing','compact');

% 1) Chord distribution
nexttile;
plot(y, span.c_y, 'LineWidth', 1.5);
grid on; xlabel('y (m)'); ylabel('c(y) (m)');
title('Chord c(y)');

% 2) Bending stiffness EI
nexttile;
plot(y, span.EI, 'LineWidth', 1.5);
grid on; xlabel('y (m)'); ylabel('EI (N·m^2)');
title('Bending stiffness EI');

% 3) Torsional stiffness GJ
nexttile;
plot(y, span.GJ, 'LineWidth', 1.5);
grid on; xlabel('y (m)'); ylabel('GJ (N·m^2)');
title('Torsional stiffness GJ');

% 4) Mass per unit length m
nexttile;
plot(y, span.m, 'LineWidth', 1.5);
grid on; xlabel('y (m)'); ylabel('m (kg/m)');
title('Mass per unit length m');

% 5) Torsion reference inertia I_sc (structural)
nexttile;
plot(y, span.I_sc, 'LineWidth', 1.5);
grid on; xlabel('y (m)'); ylabel('I_{sc} (m^4)');
title('I_{sc} (torsion ref.)');

% 6) Mass polar inertia about torsion ref (rho*I_sc)
% (you store it as span.I_sc_mass)
nexttile;
plot(y, span.I_sc_mass, 'LineWidth', 1.5);
grid on; xlabel('y (m)'); ylabel('\rho I_{sc} (kg·m)');
title('Mass inertia per length (\rho I_{sc})');

% 7) x_cm and x_sc along span
nexttile;
plot(y, span.x_cm, 'LineWidth', 1.5); hold on;
plot(y, span.x_sc, 'LineWidth', 1.5);
grid on; xlabel('y (m)'); ylabel('x (m)');
title('x_{cm}(y) and x_{sc}(y)');
legend('x_{cm}','x_{sc}','Location','best');

% 8) d = x_cm - x_sc
nexttile;
plot(y, span.d, 'LineWidth', 1.5);
grid on; xlabel('y (m)'); ylabel('d (m)');
title('Offset d = x_{cm} - x_{sc}');

% 9) dx_sc shift you applied to align SC line
%nexttile;
%plot(y, span.dx_sc, 'LineWidth', 1.5);
%grid on; xlabel('y (m)'); ylabel('\Deltax_{sc} (m)');
%title('Applied shift \Deltax_{sc}(y)');

% 9) Check shear center line perpendicular to span (x_sc_global constant)
nexttile;
x_sc_global = span.x_sc + span.dx_sc;   % global x-position of shear center line
plot(y, x_sc_global, 'LineWidth', 1.5);
grid on; xlabel('y (m)'); ylabel('x_{sc}^{glob} (m)');
title('Shear center line check: x_{sc}^{glob}(y)');

sgtitle('Cross-section properties along the span');

end