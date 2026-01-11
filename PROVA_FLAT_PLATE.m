%% flatplate_divergence_benchmark.m
% Single-file benchmark for divergence speed of a rectangular Kirchhoff-Love plate
% Requires: buildAIC.m, plateMass.m, plateStiffness.m, plateForce.m in MATLAB path

clear; clc; close all;

%% BENCHMARK INPUTS (from lecture example)
c   = 0.75;      % chord [m]
b   = 2.0;       % span [m]
h   = 0.01;      % thickness [m]
E   = 70e9;      % Young's modulus [Pa]
nu  = 0.3;       % Poisson [-]
rho = 2300;      % density [kg/m^3]

Nx = 20;         % chordwise panels
Ny = 40;         % spanwise panels

rho_inf = 1.25;  % air density [kg/m^3]
k = 0;           % reduced frequency (steady)
M_inf = 0;       % incompressible

nm = 20;         % modal subspace size for Ud (can test 10/15/20)

%% 1) STRUCTURAL MESH (exact rectangle: x in [0,c], y in [0,b])
x = linspace(0,c,Nx+1);
y = linspace(0,b,Ny+1);

N_nodes  = (Nx+1)*(Ny+1);
N_panels = Nx*Ny;
nodes = zeros(N_nodes,2);

for i = 1:Nx+1
    for j = 1:Ny+1
        kNode = (Ny+1)*(i-1) + j;
        nodes(kNode,1) = x(i);
        nodes(kNode,2) = y(j);
    end
end

% Element connectivity (same ordering as your/prof code)
elem = zeros(N_panels,4);
for i = 1:Nx
    for j = 1:Ny
        e = Ny*(i-1) + j;
        elem(e,1) = (Ny+1)*(i-1) + j;
        elem(e,2) = (Ny+1)*i + j;
        elem(e,3) = (Ny+1)*i + j+1;
        elem(e,4) = (Ny+1)*(i-1) + j+1;
    end
end

%% 2) AERODYNAMIC MESH (quarter-chord endpoints + 3/4-chord collocation)
x_p = zeros(N_panels,3);
y_p = zeros(N_panels,3);

for e = 1:N_panels
    n1 = elem(e,1); n2 = elem(e,2); n3 = elem(e,3); n4 = elem(e,4);

    x_p(e,1) = nodes(n1,1) + (nodes(n2,1) - nodes(n1,1))/4;
    y_p(e,1) = nodes(n1,2);

    x_p(e,2) = nodes(n4,1) + (nodes(n3,1) - nodes(n4,1))/4;
    y_p(e,2) = nodes(n4,2);

    x_p(e,3) = nodes(n1,1) + 3*(nodes(n2,1) - nodes(n1,1))/4;
    y_p(e,3) = 0.5*(nodes(n1,2) + nodes(n4,2));
end

%% 3) AIC (self-contained using w_doublet) - HALF SPAN + IMAGE
AIC = zeros(N_panels,N_panels);

for i = 1:N_panels
    x_i = x_p(i,3);
    y_i = y_p(i,3);

    for j = 1:N_panels
        AIC(i,j) = AIC(i,j) + w_doublet(x_i, y_i, x_p(j,:),       y_p(j,:),        c/2, M_inf, k, 1);
        AIC(i,j) = AIC(i,j) + w_doublet(x_i, y_i, x_p(j,[2,1,3]), -y_p(j,[2,1,3]),  c/2, M_inf, k, 1);
    end
end

%% 4) ASSEMBLE STRUCTURAL MATRICES (K, M, S, Ix)
N_dof = 3*N_nodes;
Mmat = zeros(N_dof,N_dof);
Kmat = zeros(N_dof,N_dof);
Smat = zeros(N_dof,N_panels);
Ix   = zeros(N_panels,N_dof);

for e = 1:N_panels
    % Element half-sizes
    a_e = 0.5*(nodes(elem(e,2),1) - nodes(elem(e,1),1));
    b_e = 0.5*(nodes(elem(e,4),2) - nodes(elem(e,1),2));

    % Element matrices
    M_e = plateMass(a_e,b_e,h,rho);
    K_e = plateStiffness(a_e,b_e,h,E,nu);
    S_e = plateForce(a_e,b_e);

    % Ix interpolation at collocation (as in lecture template)
    Ix_e = [ -9/(32*a_e), -(9*b_e)/(64*a_e),  5/32, ...
              9/(32*a_e),  (9*b_e)/(64*a_e), -3/32, ...
              9/(32*a_e), -(9*b_e)/(64*a_e), -3/32, ...
             -9/(32*a_e),  (9*b_e)/(64*a_e),  5/32 ];

    % Global DOF indices
    I_dof = [
        3*(elem(e,1)-1) + (1:3)';
        3*(elem(e,2)-1) + (1:3)';
        3*(elem(e,3)-1) + (1:3)';
        3*(elem(e,4)-1) + (1:3)';
    ];

    % Assembly
    Mmat(I_dof,I_dof) = Mmat(I_dof,I_dof) + M_e;
    Kmat(I_dof,I_dof) = Kmat(I_dof,I_dof) + K_e;
    Smat(I_dof,e)     = Smat(I_dof,e) + S_e;
    Ix(e,I_dof)       = Ix(e,I_dof) + Ix_e;
end

%% 5) BOUNDARY CONDITIONS (clamped root y=0 => all DOFs at j=1 fixed)
I_fix = [ ...
    1:3*(Ny+1):N_dof, ...
    2:3*(Ny+1):N_dof, ...
    3:3*(Ny+1):N_dof ...
];
I_free = setdiff(1:N_dof, I_fix);

K_free = Kmat(I_free,I_free);
M_free = Mmat(I_free,I_free);

%% 6) A0 operator and reduced divergence solve
A0 = 0.5*rho_inf * (Smat * (AIC \ Ix));
A0_free = A0(I_free,I_free);

% modal basis (sorted + M-orthonormal)
[Phi, D] = eig(K_free, M_free);
w2 = real(diag(D));
[~,idx] = sort(w2,'ascend');
Phi = Phi(:,idx);

for j = 1:size(Phi,2)
    Phi(:,j) = Phi(:,j) / sqrt(Phi(:,j)'*M_free*Phi(:,j));
end

Phi_r = Phi(:,1:nm);
K_red  = Phi_r' * K_free  * Phi_r;
A0_red = Phi_r' * A0_free * Phi_r;

% symmetrize in reduced space for steady divergence
A0_red = 0.5*(A0_red + A0_red.');

lam = eig(K_red, A0_red);
mask = (real(lam)>0) & (abs(imag(lam))<1e-8);
Ud = sqrt(min(real(lam(mask))));

%% 7) PRINT CHECKS
fprintf('\n--- Flat plate divergence benchmark ---\n');
fprintf('c=%.3f b=%.3f h=%.3f E=%.2e nu=%.2f rho=%.0f Nx=%d Ny=%d\n', c,b,h,E,nu,rho,Nx,Ny);
fprintf('nm=%d\n', nm);
fprintf('Ud = %.4f m/s\n', Ud);

% Quick diagnostics
symA = norm(A0_free-A0_free','fro')/norm(A0_free,'fro');
fprintf('sym(A0_free) = %.4f\n', symA);
fprintf('rcond(AIC)   = %.3e\n', rcond(AIC));
