clear
close all

% Problem data
c = 0.75;
b = 2;
h = 0.01;
E = 70e9;
G = E/(2*(1+0.3));
rho = 2300;
N = 100;

%Section properties
A = h*c;
mu = rho*A;
x_cm = 0.5*c;
x_sc = 0.5*c;
I = c*h^3/12;
I_sc = rho*c^3*h/12; %At the shear center
k_t = 0.000705; % Hi ha un livescript per a calcular les section properties (com aquest param)
J = k_t*c^3*h/12;


% Structural mesh

y_ = (0:b/N:b)'; %wtf ha fet aqui


% Aerodynamic mesh

x_p = zeros(N,5);
y_p = zeros(N,5);

for i = 1:N
    x_p(i,[1,4]) = 20*c; %El punt 1 es el més llunyà. 20 cordes ha de ser suficient.
    x_p(i,[2,3]) = c/4;
    x_p(i,5) = 3*c/4;

    y_p(i,[1,2]) = y_(i);
    y_p(i,[3,4]) = y_(i+1);
    y_p(i,5) = (y_(i)+y_(i+1))/2;


end


%Aerodynamic influence coefficients

%AIC = horseshoeAIC(x_p,y_p); %D'aquesta funcio tinc foto al mobil. No tinc molt clar si s'hauria de poder obtenir del horseshoe_vortex live script

%Assembly of structural matrices

N_dof = 3*(N+1);
M = zeros(N_dof,N_dof);
K = zeros(N_dof,N_dof);
B = zeros(N_dof,N);

for i = 1:N
    % Element size
    b_e = y_(i+1)-y_(i);

    %Element mass matrix
    M_e = beamMass(b_e,mu,I_sc,x_cm,x_sc);
    %Element mass matrix
    K_e = beamStiffness(b_e,mu,I_sc,x_cm,x_sc);    
    %Element mass matrix
    B_e = beamForce(b_e,c/4,x_sc); 

    %Aquestes són les funcions que genera el codi del profe, i que fem
    %servir aqui.

    %Global indices
    I_dof = 3*([i,i,i,i+1,i+1,i+1]-1) + [1,2,3,1,2,3];

    %Assembly in global matrices
    M(I_dof,I_dof) = M(I_dof,I_dof) + M_e;
    K(I_dof,I_dof) = K(I_dof,I_dof) + K_e;
    B(I_dof,i) = B(I_dof,i) + B_e;
end


%Boundary conditions
I_free = 4:N_dof;
M_free = M(I_free,I_free);
K_free = M(I_free,I_free);
B_free = M(I_free,:);

% Structural modes
q_mod = zeros(N_dof,N_dof-3);
[q_mod(I_free,:),w2] = eig(K_free,M_free);
freq = sqrt(diag(w2))/(2*pi);

%Plot FOTO AL MOBIL




%Constant angle of attack

%Flow properties
rho_inf = 1.2;
U_inf = 20;
alpha = 15*pi/180;

%Vertical velocity vector
W_ref = U_inf*alpha*ones(N,1);

%Lift distribution
lift = rho_inf*U_inf*(-AIC\wref);

% DOFs caused by angle of attack
q_aoa = zeros(N_dof,1);
q_aoa(I_free) = K_free\(B_free*lift);


%plots ??