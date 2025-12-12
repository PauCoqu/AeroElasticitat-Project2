function res = w_doublet(x,y,x_panel,y_panel,b,M,k,Cp)
% Compute the kernels
rx_0 = (x-(x_panel(1)+x_panel(2))/2)/b;
ry_0 = (y-(y_panel(1)+y_panel(2))/2)/b;
Kbar_0 = Kbar(rx_0,ry_0,M,k);
rx_1 = (x-x_panel(1))/b;
ry_1 = (y-y_panel(1))/b;
Kbar_1 = Kbar(rx_1,ry_1,M,k);
rx_2 = (x-x_panel(2))/b;
ry_2 = (y-y_panel(2))/b;
Kbar_2 = Kbar(rx_2,ry_2,M,k);

% Compute parameters
c_k = (2*x_panel(3)-x_panel(1)-x_panel(2))/b;
b_k = sqrt((x_panel(2)-x_panel(1))^2 + (y_panel(2)-y_panel(1))^2)/b;
sinL_k = (y_panel(2)-y_panel(1))/b/b_k;

% Compute A coefficients
A_0 = Kbar_0;
A_1 = (Kbar_2-Kbar_1)/b_k;
A_2 = (2*Kbar_1-4*Kbar_0+2*Kbar_2)/b_k^2;

% Compute B coefficients
fact0 = 4*c_k*b_k/(4*ry_0^2-b_k^2*sinL_k^2);
fact1 = log((b_k*sinL_k-2*ry_0)^2/(b_k*sinL_k+2*ry_0)^2);
B_0 = fact0*A_0;
B_1 = (c_k/(2*sinL_k^2))*fact1*A_1 ...
    + (ry_0/sinL_k)*fact0*A_1;
B_2 = (c_k*b_k/sinL_k^2)*A_2 + ...
    + (ry_0*c_k/sinL_k^3)*fact1*A_2 ...
    + (ry_0^2/sinL_k^2)*fact0*A_2;

% Compute induced velocity
res = -Cp/(8*pi)*(B_0 + B_1 + B_2);
end