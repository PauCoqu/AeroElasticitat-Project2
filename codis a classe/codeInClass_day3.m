%% Parts reused from other codes
clear
close all;

c=0.75;
b=2;
h =0.01;
xh = 0.6;
y1 = 1.25;
y2= 1.75;
E = 70e9;
nu = 0.3;
rho = 2300;

Nx = 20;
Ny = 40;


% Structural two number of elements, to account for the control surface and
% three for ny
Nx2 = ceil(Nx*(c-xh)/c);
Nx1 = Nx - Nx2;

Ny3 = ceil(Ny*(b-y2)/b);
Ny2 = ceil(Ny*(y2-y1)/b);
Ny1 = Ny - Ny2 - Ny3;

x = [0:xh/Nx1:xh,xh+(c-xh)/Nx2:(c-xh)/Nx2:c]; %Crec que aixo esta mal

y = [0:y1/Ny1:y1,y1+(y2-y1)/Ny2:(y2-y1(Ny2:y2)),y2+(b-y2)/Ny3:(b-y2)/Ny3:b];

nodes = zeros((Nx+1)*(Ny+1),2);

for i = 1:Nx+1
    for j = 1:Ny+1
        k = (Ny+1)*(i-1) + j;
        nodes(k,1) = x(i);
        nodes(k,2) = y(j);
    end
end


% Element nodes
elem = zeros(Nx*Ny,4);
control = false(Nx*Ny,4);
for i = 1:Nx
    for j = 1:Ny
        e = Ny*(i-1) + j;
        elem(e,1) = (Ny+1)*(i-1) + j;
        elem(e,2) = (Ny+1)*i + j;
        elem(e,3) = (Ny+1)*i + j+1;
        elem(e,4) = (Ny+1)*(i-1) + j+1;

        if i > Nx1 && j > Ny1 &&  j <= Ny1+Ny2 % Check if panel is in control surface or not
            control(e) = true;
        end
    end
end

% Control surface nodes
i_delta = unique(elem(control,:)); % I dont know how but it should only output the true ones

%% AERODYNAMIC MESH

x_p = zeros(Nx*Ny,3);
y_p = zeros(Nx*Ny,3);

for e = 1:Nx*Ny
    x_p(e,1) = nodes(elem(e,1),1) + (nodes(elem(e,2),1) - nodes(elem(e,1),1))/4;
    y_p(e,1) = nodes(elem(e,1),2); 
    x_p(e,2) = nodes(elem(e,4),1) + (nodes(elem(e,3),1) - nodes(elem(e,4),1))/4;
    y_p(e,2) = nodes(elem(e,4),2);
    x_p(e,3) = nodes(elem(e,1),1) + 3*(nodes(elem(e,2),1) - nodes(elem(e,1),1))/4;
    y_p(e,3) = (nodes(elem(e,1),2) + nodes(elem(e,4),2))/2;
end


%AICs
k=0; % Steady case
M_inf = 0; % Assume incompressibility conditions

AIC = buildAIC(x_p,y_p,k,M_inf,c/2);


% Assembly of the structural matrices
%Assembly of structural matrices
N_dof = 3*size(nodes,1);
M = zeros(N_dof,N_dof);
K = zeros(N_dof,N_dof);
S = zeros(N_dof,Nx*Ny);
Ix = zeros(Nx*Ny,N_dof);

for i = 1:Nx*Ny
    % Element size
    a_e = (nodes(elem(i,2),1) - nodes(elem(i,1),1))/2;
    b_e = (nodes(elem(i,4),2) - nodes(elem(i,1),2))/2;
    % Element matrices
    M_e = plateMass(a_e,b_e,h,rho);
    K_e = plateStiffness(a_e,b_e,h,E,nu);
    S_e = plateForce(a_e,b_e);
    % Interpolation matrix collocation point
    Ix_e = [
            -9/(32*a_e), -9*b_e/64*a_e, 5/32,
            9/(32*a_e), 9*b_e/64*a_e, -5/32,
            9/(32*a_e), -9*b_e/64*a_e, 5/32,
            -9/(32*a_e), 9*b_e/64*a_e, -5/32,

        ]


    %Globsl indices vector
    I_dof = [
        3*(elem(i,1)-1) + [1;2;3];
        3*(elem(i,2)-1) + [1;2;3];
        3*(elem(i,3)-1) + [1;2;3];
        3*(elem(i,4)-1) + [1;2;3];
        ];

    %Assembly
    M(I_dof,I_dof) = M(I_dof,I_dof) + M_e;
    K(I_dof,I_dof) = K(I_dof,I_dof) + K_e;
    S(I_dof,i) = S(I_dof,i) +  S_e;
    Ix(i,I_dof) = Ix(i,I_dof) + Ix_e;
end

%Boundary conditions
I_fix = [
    1:3*(Ny+1):N_dof,...
    2:3*(Ny+1):N_dof,...
    3:3*(Ny+1):N_dof
    ];
I_free = setdiff(1:N_dof,I_fix);

%Reduced matrices
M_free = M(I_free,I_free);
K_free = K(I_free,I_free);
S_free = S(I_free,:);

% Aeroelastic coupling
rho_inf = 1.25;
A0 = rho_inf/2*S*(AIC\Ix);

% Comput divergence speed
[X_div,U_inf] = eigs(sparse(K(I_free,I_free)),sparse(A0(I_free,I_free)),10,"sm");



%Structural
q_mod = zeros(N_dof,length(I_free));

[q_mod(I_free,:),w2] = eig(K_free,M_free);
freq = sqrt(diag(w2))/(2*pi);

% Reduced system
i_modes = 1:4;
K_red = q_mod(:,i_modes)'*K*q_mod(:,i_modes);
A0_red = q_mod(:,i_modes)'*A0*q_mod(:,i_modes);

[X_div,U2_inf] = eig(K_red,A0_red);
U_inf = sort(diag(sqrt(U2_inf)));
U_d = U_inf(1);


% FOTO DE CLASSE AL MBOIL

% Loop through the velocities


for i = 1:lenght(U_inf)
    % Eff stiffness matrix
    K_eff = K_red - U_inf(i)^2*A0_red;

    % Plate deformation
    q(:,i) = q_mod(:,i_modes)*(K_eff\(U_inf(i)^2*q_mod(:,i_modes)'*A0*q_delta));

    % Pressure distribution
    delta_p(:,i) = 0.5*rho_inf*U_inf(i)^2*(AIC\Ix*(q(:,i) + q_delta)));
    delta_p_0(:,i) = 0.5*rho_inf*U_inf(i)^2*(AIC\Ix*q_delta));
end


% Compute the total lift
L = sum(S(1:3:end,:)*delta_p,1);
L_0 = sum(S(1:3:end,:)*delta_p_0,1);


%Plot




