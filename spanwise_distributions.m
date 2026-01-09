function [span] = spanwise_distributions(section, material, c_y, c_root, y, b, S, AR, lambda)

A0   = section.A;
Ixx0 = section.Ixx;
x_cm0  = section.xcm;
x_sc0  = section.xsc;      
I_sc0  = section.Isc;


% Scaling 
A   = A0  * (c_y/c_root).^2;
Ixx = Ixx0* (c_y/c_root).^4;
I_sc = I_sc0 * (c_y/c_root).^4;
J   = section.kt * I_sc;
x_cm = x_cm0(1)*(c_y/c_root);
x_sc = x_sc0(1)*(c_y/c_root);
I_sc_mass = material.Density.*I_sc;

% Material properties
E  = material.YoungModulus;
nu = material.Poisson;
G  = E/(2*(1+nu));
EI = E * Ixx;
GJ = G * J;
m  = material.Density * A;

% Distance between CM and SC
d  = x_cm - x_sc;

% Distance moved
dx_sc = x_sc0(1)-x_sc; % cant forget the distance we have moved it

% Store everything
span.y      = y;
span.c_y    = c_y;
span.EI     = EI;
span.GJ     = GJ;
span.m      = m;
span.I_sc   = I_sc;
span.I_sc_mass = I_sc_mass;
span.x_cm   = x_cm;
span.x_sc   = x_sc;
span.d      = d;
span.b      = b;
span.S      = S;
span.c_root = c_root;
span.lambda = lambda;
span.AR     = AR;
span.dx_sc  = dx_sc; % needed distance to obtain a shear center line perpendicular to plane xz.

end
