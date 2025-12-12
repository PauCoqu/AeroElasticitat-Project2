function [M_e, K_e, B_e] = beamElementSymbolic()
%% 3.1 Bending
syms t eta_hat_1(t) eta_hat_2(t) zeta_hat_1(t) zeta_hat_2(t)

% Let's define the coordinates mapping
syms xi y_hat_1 y_hat_2 b_e
y = (y_hat_1+y_hat_2)/2 + b_e/2*xi;  %ABANS DEIA: y(xi) = (y_hat_1+y_hat_2)/2 + b_e/2*xi;  

% Cubic Hermite shape functions definition in terms of xi
N_eta_1(xi) = (2-3*xi+xi^3)/4;
N_zeta_1(xi) = (1-xi-xi^2+xi^3)/4;
N_eta_2(xi) = (2+3*xi-xi^3)/4;
N_zeta_2(xi) = (-1-xi+xi^2+xi^3)/4;
figure("Units","normalized","Position",[0,0,1,1])
subplot(2,2,1)
fplot(N_eta_1,[-1,1]); axis equal; ylim([0,1]); xlabel \xi; title N_{1}^{[\eta]};
subplot(2,2,3)
fplot(N_zeta_1,[-1,1]); axis equal; ylim([0,1]); xlabel \xi; title N_{1}^{[\zeta]};
subplot(2,2,2)
fplot(N_eta_2,[-1,1]); axis equal; ylim([0,1]); xlabel \xi; title N_{2}^{[\eta]};
subplot(2,2,4)
fplot(N_zeta_2,[-1,1]); axis equal; ylim([-1,0]); xlabel \xi; title N_{2}^{[\zeta]};

% Deflection expression inside the element
eta(xi,t) = N_eta_1(xi)*eta_hat_1(t) + N_zeta_1(xi)*b_e/2*zeta_hat_1(t) ...
          + N_eta_2(xi)*eta_hat_2(t) + N_zeta_2(xi)*b_e/2*zeta_hat_2(t);

% The derivatives of the deflections can be directly obtained
zeta(xi,t) = 2/b_e*diff(eta(xi,t),xi);
dzeta_dy(xi,t) = 2/b_e*diff(zeta(xi,t),xi);
deta_dt(xi,t) = diff(eta(xi,t),t);

%% 3.2 Twist

% Let's define the twist DOFs and shape functions
syms theta_hat_1(t) theta_hat_2(t)
N_theta_1(xi) = (1-xi)/2;
N_theta_2(xi) = (1+xi)/2;
figure("Units","normalized","Position",[0,0,1,0.5])
subplot(1,2,1)
fplot(N_theta_1,[-1,1]); axis equal; ylim([0,1]); xlabel \xi; title N_{1}^{[\theta]};
subplot(1,2,2)
fplot(N_theta_2,[-1,1]); axis equal; ylim([0,1]); xlabel \xi; title N_{2}^{[\theta]};
% Twist along the element
theta(xi,t) = N_theta_1(xi)*theta_hat_1(t) + N_theta_2(xi)*theta_hat_2(t);

% The derivatives of the twist can be directly obtained
dtheta_dy(xi,t) = 2/b_e*diff(theta(xi,t),xi);
dtheta_dt(xi,t) = diff(theta(xi,t),t);

%% 3.3 Energy functions

% Definition of the mass properties of the cross-section
syms mu_e I_sc_e x_cm_e x_sc
% Definition of the elastic properties of the cross-section
syms E_e I_e G_e J_e
% Definition of the aerodynamic properties
syms l_e x_ac_e

% Kinetic energy definition
T_e = mu_e*b_e/2*int(deta_dt(xi,t)^2/2,xi,-1,1) ...
    + I_sc_e*b_e/2*int(dtheta_dt(xi,t)^2/2,xi,-1,1) ...
    + mu_e*b_e/2*(x_sc-x_cm_e)*int(deta_dt(xi,t)*dtheta_dt(xi,t),xi,-1,1);

% Potential energy definition
U_e = E_e*I_e*b_e/2*int(dzeta_dy(xi,t)^2/2,xi,-1,1) ...
    + G_e*J_e*b_e/2*int(dtheta_dy(xi,t)^2/2,xi,-1,1);

% Virtual work definition
dW_e = l_e*b_e/2*int(eta(xi,t) + (x_sc-x_ac_e)*theta(xi,t),xi,-1,1);

%% 3.4 Determination of the element matrices

% Let's define our DOFs vector as a cell array (since the components are symbolic functions, it is better to define it as a cell array)
q = {eta_hat_1, zeta_hat_1, theta_hat_1, eta_hat_2, zeta_hat_2, theta_hat_2};

% Matrices definition
M_e = sym('M_e',[6,6]);
K_e = sym('K_e',[6,6]);
B_e = sym('B_e',[6,1]);
for i = 1:6
    for j = 1:6
        M_e(i,j) = diff(diff(diff(T_e,diff(q{i},t)),t),diff(q{j},t,2));
        K_e(i,j) = diff(diff(U_e,q{i}),q{j});
    end
    B_e(i) = diff(diff(dW_e,q{i}),l_e);
end

% Now we can transform the element matrices symbolic expressions into
% MATLAB functions that can then be used to generate these matrices

matlabFunction(M_e,'File','beamMass.m','vars',[b_e,mu_e,I_sc_e,x_cm_e,x_sc],'Optimize',true);
matlabFunction(K_e,'File','beamStiffness.m','vars',[b_e,E_e,I_e,G_e,J_e],'Optimize',true);
matlabFunction(B_e,'File','beamForce.m','vars',[b_e,x_ac_e,x_sc],'Optimize',true);

end