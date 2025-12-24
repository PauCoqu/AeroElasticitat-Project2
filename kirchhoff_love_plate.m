
% We define the DOFs at each node as symbolic variables
syms t eta_hat_1(t) eta_hat_2(t) eta_hat_3(t) eta_hat_4(t) ...
       zeta_hat_1(t) zeta_hat_2(t) zeta_hat_3(t) zeta_hat_4(t) ...
       theta_hat_1(t) theta_hat_2(t) theta_hat_3(t) theta_hat_4(t)

% Cubic Hermite shape functions definition in terms of xi
syms xi
H_eta_1(xi) = (2-3*xi+xi^3)/4;
H_zeta_1(xi) = (1-xi-xi^2+xi^3)/4;
H_eta_2(xi) = (2+3*xi-xi^3)/4;
H_zeta_2(xi) = (-1-xi+xi^2+xi^3)/4;

% Bi-cubic Hermite shape functions
syms xi_1 xi_2
N_eta_1(xi_1,xi_2) = H_eta_1(xi_1)*H_eta_1(xi_2);
N_zeta_1(xi_1,xi_2) = H_eta_1(xi_1)*H_zeta_1(xi_2);
N_theta_1(xi_1,xi_2) = -H_zeta_1(xi_1)*H_eta_1(xi_2);
N_eta_2(xi_1,xi_2) = H_eta_2(xi_1)*H_eta_1(xi_2);
N_zeta_2(xi_1,xi_2) = H_eta_2(xi_1)*H_zeta_1(xi_2);
N_theta_2(xi_1,xi_2) = -H_zeta_2(xi_1)*H_eta_1(xi_2);
N_eta_3(xi_1,xi_2) = H_eta_2(xi_1)*H_eta_2(xi_2);
N_zeta_3(xi_1,xi_2) = H_eta_2(xi_1)*H_zeta_2(xi_2);
N_theta_3(xi_1,xi_2) = -H_zeta_2(xi_1)*H_eta_2(xi_2);
N_eta_4(xi_1,xi_2) = H_eta_1(xi_1)*H_eta_2(xi_2);
N_zeta_4(xi_1,xi_2) = H_eta_1(xi_1)*H_zeta_2(xi_2);
N_theta_4(xi_1,xi_2) = -H_zeta_1(xi_1)*H_eta_2(xi_2);

% Plot shape functions
figure("Units","normalized","Position",[0,0,0.5,1])
subplot(4,3,1)
fsurf(N_eta_1,[-1,1,-1,1],"EdgeColor","none"); axis equal; zlim([0,1]);
xlabel \xi_1; ylabel \xi_2; title N_{1}^{[\eta]};
subplot(4,3,2)
fsurf(N_zeta_1,[-1,1,-1,1],"EdgeColor","none"); axis equal; zlim([0,1]);
xlabel \xi_1; ylabel \xi_2; title N_{1}^{[\zeta]};
subplot(4,3,3)
fsurf(N_theta_1,[-1,1,-1,1],"EdgeColor","none"); axis equal; zlim([-1,0]);
xlabel \xi_1; ylabel \xi_2; title N_{1}^{[\theta]};
subplot(4,3,4)
fsurf(N_eta_2,[-1,1,-1,1],"EdgeColor","none"); axis equal; zlim([0,1]);
xlabel \xi_1; ylabel \xi_2; title N_{2}^{[\eta]};
subplot(4,3,5)
fsurf(N_zeta_2,[-1,1,-1,1],"EdgeColor","none"); axis equal; zlim([0,1]);
xlabel \xi_1; ylabel \xi_2; title N_{2}^{[\zeta]};
subplot(4,3,6)
fsurf(N_theta_2,[-1,1,-1,1],"EdgeColor","none"); axis equal; zlim([0,1]);
xlabel \xi_1; ylabel \xi_2; title N_{2}^{[\theta]};
subplot(4,3,7)
fsurf(N_eta_3,[-1,1,-1,1],"EdgeColor","none"); axis equal; zlim([0,1]);
xlabel \xi_1; ylabel \xi_2; title N_{3}^{[\eta]};
subplot(4,3,8)
fsurf(N_zeta_3,[-1,1,-1,1],"EdgeColor","none"); axis equal; zlim([-1,0]);
xlabel \xi_1; ylabel \xi_2; title N_{3}^{[\zeta]};
subplot(4,3,9)
fsurf(N_theta_3,[-1,1,-1,1],"EdgeColor","none"); axis equal; zlim([0,1]);
xlabel \xi_1; ylabel \xi_2; title N_{3}^{[\theta]};
subplot(4,3,10)
fsurf(N_eta_4,[-1,1,-1,1],"EdgeColor","none"); axis equal; zlim([0,1]);
xlabel \xi_1; ylabel \xi_2; title N_{4}^{[\eta]};
subplot(4,3,11)
fsurf(N_zeta_4,[-1,1,-1,1],"EdgeColor","none"); axis equal; zlim([-1,0]);
xlabel \xi_1; ylabel \xi_2; title N_{4}^{[\zeta]};
subplot(4,3,12)
fsurf(N_theta_4,[-1,1,-1,1],"EdgeColor","none"); axis equal; zlim([-1,0]);
xlabel \xi_1; ylabel \xi_2; title N_{4}^{[\theta]};

% Then, we can define the deflection and its derivatives:

% Deflection expression inside the element
syms a_e b_e
eta(xi_1,xi_2,t) = ...
    N_eta_1(xi_1,xi_2)*eta_hat_1(t) + ...
    N_zeta_1(xi_1,xi_2)*b_e*zeta_hat_1(t) + ...
    N_theta_1(xi_1,xi_2)*a_e*theta_hat_1(t) + ...
    N_eta_2(xi_1,xi_2)*eta_hat_2(t) + ...
    N_zeta_2(xi_1,xi_2)*b_e*zeta_hat_2(t) + ...
    N_theta_2(xi_1,xi_2)*a_e*theta_hat_2(t) + ...
    N_eta_3(xi_1,xi_2)*eta_hat_3(t) + ...
    N_zeta_3(xi_1,xi_2)*b_e*zeta_hat_3(t) + ...
    N_theta_3(xi_1,xi_2)*a_e*theta_hat_3(t) + ...
    N_eta_4(xi_1,xi_2)*eta_hat_4(t) + ...
    N_zeta_4(xi_1,xi_2)*b_e*zeta_hat_4(t) + ...
    N_theta_4(xi_1,xi_2)*a_e*theta_hat_4(t);

% Derivatives
deta_dx(xi_1,xi_2,t) = 1/a_e*diff(eta(xi_1,xi_2,t),xi_1);
deta_dy(xi_1,xi_2,t) = 1/b_e*diff(eta(xi_1,xi_2,t),xi_2);
d2eta_dx2(xi_1,xi_2,t) = 1/a_e*diff(deta_dx(xi_1,xi_2,t),xi_1);
d2eta_dy2(xi_1,xi_2,t) = 1/b_e*diff(deta_dy(xi_1,xi_2,t),xi_2);
d2eta_dxdy(xi_1,xi_2,t) = 1/(a_e*b_e)*diff(diff(eta(xi_1,xi_2,t),xi_1),xi_2);
deta_dt(xi_1,xi_2,t) = diff(eta(xi_1,xi_2,t),t);

%% 2. Energy functions
% 2.1. Kinetic energy

% Definition of the mass properties
syms rho_e h_e

% Kinetic energy definition
T_e = rho_e*h_e*a_e*b_e*int(int(deta_dt(xi_1,xi_2,t)^2/2,xi_1,-1,1),xi_2,-1,1);

% 2.2. Potential energy
% Definition of the stiffness properties
syms E_e nu_e

% Potential energy definition
U_e = E_e*h_e^3*a_e*b_e/(24*(1-nu_e^2))*int(int(...
    d2eta_dx2(xi_1,xi_2,t)^2 + d2eta_dy2(xi_1,xi_2,t)^2 + ...
    2*nu_e*d2eta_dx2(xi_1,xi_2,t)*d2eta_dy2(xi_1,xi_2,t) + ...
    2*(1-nu_e)*d2eta_dxdy(xi_1,xi_2,t)^2,xi_1,-1,1),xi_2,-1,1);

% 2.3. Virtual work
% Definition of the aerodynamic properties
syms Delta_p_e

% Virtual work definition
dW_e = Delta_p_e*a_e*b_e*int(int(eta(xi_1,xi_2,t),xi_1,-1,1),xi_2,-1,1);

%% 3. Determination of the element matrices

% Let's define our DOFs vector as a cell array (since the components are
% symbolic functions, it is better to define it as a cell array)
q = {eta_hat_1, zeta_hat_1, theta_hat_1, ...
     eta_hat_2, zeta_hat_2, theta_hat_2, ...
     eta_hat_3, zeta_hat_3, theta_hat_3, ...
     eta_hat_4, zeta_hat_4, theta_hat_4};

% Matrices definition
M_e = sym('M_e',[12,12]);
K_e = sym('K_e',[12,12]);
S_e = sym('S_e',[12,1]);
for i = 1:12
    for j = 1:12
        M_e(i,j) = diff(diff(diff(T_e,diff(q{i},t)),t),diff(q{j},t,2));
        K_e(i,j) = diff(diff(U_e,q{i}),q{j});
    end
    S_e(i) = diff(diff(dW_e,q{i}),Delta_p_e);
end

% Now we can transform the element matrices symbolic expressions into
% MATLAB functions that can then be used to generate these matrices
matlabFunction(M_e,'File','plateMass.m',...
                   'vars',[a_e,b_e,h_e,rho_e],'Optimize',true);
matlabFunction(K_e,'File','plateStiffness.m',...
                   'vars',[a_e,b_e,h_e,E_e,nu_e],'Optimize',true);
matlabFunction(S_e,'File','plateForce.m',...
                   'vars',[a_e,b_e],'Optimize',true);