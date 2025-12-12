function [span] = spanwise_distributions(section, material, c_y, c_root, y, b, S, AR, lambda)

A0   = section.A;
Ixx0 = section.Ixx;
J0   = section.kt;
x_cm0  = section.xcm;
x_sc0  = section.xsc;      
I_sc0  = section.Isc;      

% Scaling 
A   = A0  * (c_y/c_root).^2;
Ixx = Ixx0* (c_y/c_root).^4;
J   = J0  * (c_y/c_root).^4;
I_sc = I_sc0 * (c_y/c_root).^4;          
x_cm = x_cm0(1)*(c_y/c_root);
x_sc = x_sc0(1)*(c_y/c_root);

% Material properties
E  = material.YoungModulus;
nu = material.Poisson;
G  = E/(2*(1+nu));
EI = E * Ixx;
GJ = G * J;
m  = material.Density * A;

% Distance between CM and SC
d  = x_cm - x_sc;

% Store everything
span.y      = y;
span.c_y    = c_y;
span.EI     = EI;
span.GJ     = GJ;
span.m      = m;
span.I_sc   = I_sc;
span.x_cm   = x_cm;
span.x_sc   = x_sc;
span.d      = d;
span.b      = b;
span.S      = S;
span.c_root = c_root;
span.lambda = lambda;
span.AR     = AR;
end
