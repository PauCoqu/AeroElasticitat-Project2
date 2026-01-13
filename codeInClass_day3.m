%% ES L'EXEMPLE DE CLASSE DE UD = 145 APROX (FLAT PLATE)

clear; clc; close all

c=1;
b=2;
h =0.01;
xh = 0.6;
y1 = 1.25;
y2= 1.5;
E = 70e9;
nu = 0.3;
rho = 2300;
rho_inf = 1.25;

Nx = 20;
Ny = 40;


% Structural two number of elements, to account for the control surface and
% three for ny
Nx2 = ceil(Nx*(c-xh)/c);
Nx1 = Nx - Nx2;

Ny3 = ceil(Ny*(b-y2)/b);
Ny2 = ceil(Ny*(y2-y1)/b);
Ny1 = Ny - Ny2 - Ny3;

x = [0:xh/Nx1:xh,xh+(c-xh)/Nx2:(c-xh)/Nx2:c];
y = [0:y1/Ny1:y1,y1+(y2-y1)/Ny2:(y2-y1)/Ny2:y2,y2+(b-y2)/Ny3:(b-y2)/Ny3:b];

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
It = zeros(Nx*Ny,N_dof);

for i = 1:Nx*Ny
    % Element size
    a_e = (nodes(elem(i,2),1) - nodes(elem(i,1),1))/2;
    b_e = (nodes(elem(i,4),2) - nodes(elem(i,1),2))/2;
    % Element matrices
    M_e = plateMass(a_e,b_e,h,rho);
    K_e = plateStiffness(a_e,b_e,h,E,nu);
    S_e = plateForce(a_e,b_e);
    %Coeficients del professor (el profe tenia 5/32 al 2 i al 3, teoria posa 3/32)
   Ix_e = [ -9/(32*a_e), -(9*b_e)/(64*a_e),  5/32, ...
             9/(32*a_e),  (9*b_e)/(64*a_e), -3/32, ...
             9/(32*a_e), -(9*b_e)/(64*a_e), -3/32, ...
            -9/(32*a_e),  (9*b_e)/(64*a_e),  5/32 ];

   
    % Collocation point interpolation
    It_e = [5/64, 5*b_e/128, -3*a_e/64, ...
            27/64, 27*b_e/128, 9*a_e/64, ...
            27/64, -27*b_e/128, 9*a_e/64, ...
            5/64, -5*b_e/128, -3*a_e/64];


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
    It(i,I_dof) = It(i,I_dof) + It_e;
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
Ix_free = Ix(:,I_free);
It_free = It(:,I_free);

% Aeroelastic coupling
A0 = rho_inf/2*S*(AIC\Ix);

% Comput divergence speed ALL MODES
[X_div,U2_inf] = eigs(sparse(K(I_free,I_free)), sparse(A0(I_free,I_free)), 10, 'smallestabs');
U_inf_ALLMODES = sort(diag(sqrt(U2_inf)));
Ud_allmodes = U_inf_ALLMODES(1) %DIVERGENCE SPEED

%Structural
q_mod = zeros(N_dof,length(I_free));
[q_mod(I_free,:),w2] = eig(K_free,M_free);
freq = sqrt(diag(w2))/(2*pi); %frequencia

% Reduced system 
i_modes = 1:3; %Si ometem el primer mode baixa molt (144...)
K_red = q_mod(:,i_modes)'*K*q_mod(:,i_modes);
A0_red = q_mod(:,i_modes)'*A0*q_mod(:,i_modes);

[X_div,U2_inf2] = eig(K_red,A0_red);
U_inf_REDUCED = sort(diag(sqrt(U2_inf2)));
Ud_reduced = U_inf_REDUCED(1) %DIVERGENCE SPEED

%% AFEGIT A PARTIR D'AQUÍ
%% 2. RISK OF FLUTTER (P METHOD)
%k = 0; % quasi-steady  %JA DEFINIT MÉS AMUNT
%M_inf = 0; % Incompressibility %JA DEFINIT MÉS AMUNT
Ud = Ud_allmodes; %Aprop de Ud = 100 ens donen uns grafics prou correctes

% Select modes for p-method as we do a Model reduction
c_root = c;
i_modes = [1,2];
N_modes = length(i_modes);
% Previous reduced matrices
% Phi = q_mod(:, i_modes); PREVIOUS ONE USED
% M_red  = Phi.'*M*Phi; PREVIOUS ONE USED
% K_red  = Phi.'*K*Phi; PREVIOUS ONE USED
% S_red  = Phi.'*S; PREVIOUS ONE USED
% Ix_red = Ix*Phi; PREVIOUS ONE USED
% It_red = It*Phi; PREVIOUS ONE USED
% A0_red = rho_inf/2*S_red*(AIC\Ix_red); PREVIOUS ONE USED
% A1_red = -rho_inf/2*S_red*(AIC\It_red); PREVIOUS ONE USED


% Reduced matrices
Phi = q_mod(I_free, i_modes);   %use free part of the modes WHY?
M_red  = Phi.'*M_free*Phi;
K_red  = Phi.'*K_free*Phi;
A0_free = rho_inf/2 * S_free * (AIC \ Ix_free);
A1_free = -rho_inf/2 * S_free * (AIC \ It_free);
A0_red = Phi.'*A0_free*Phi;
A1_red = Phi.'*A1_free*Phi;

% Initialize velocity vector
dU = 5;
U_ = 5:5:1.6*Ud;

% Zero tolerance
tol = 1e-6;

% Initialize variables
p_ = nan(2*N_modes,length(U_));
m_ = nan(N_modes,2*N_modes,length(U_));
U_min = [];
U_max = [];

% Loop through velocities
for i = 1:length(U_)

    % Effecftive matrices
    Keff = K_red - U_(i)^2*A0_red;
    Ceff = U_(i)*A1_red;
    Meff = M_red;


    % Generalized eigenvalues
    O = zeros(size(Keff));
    I = eye(size(Keff));
    A = [Keff, O; O, I];
    B = [-Ceff, -Meff; I, O];
    [X,P] = eig(A,B);
    p = diag(P);

    % Sorting algorithm
    if i == 1
        p_(:,i) = p;
        m_(:,:,i) = X(1:N_modes,:);
    else
        m2sort = 1:2*N_modes;
        for j = 1:length(p)
            [~,jmin] = min(abs(real(p_(j,i-1))-real(p(m2sort)))+ abs(imag(p_(j,i-1))-imag(p(m2sort))));

            p_(j,i) = p(m2sort(jmin));
            m_(:,j,i) = X(1:N_modes,m2sort(jmin));
            m2sort = setdiff(m2sort,m2sort(jmin));

        end
    end

if i > 1
    r_prev = max(real(p_(:,i-1)));
    if isempty(U_min) && (r_prev <= 0) && (max(real(p_(:,i))) > 0)
        U_min = U_(i-1)
        U_max = U_(i)
    end
end
end

% PLOTS

figure
subplot(2,1,1)
hold on; box on;
plot(U_/Ud, real(p_).*c_root./(2*U_), 'LineWidth', 2);
xlabel("U_{\infty}/U_D");
ylabel("p_R c / 2U");

subplot(2,1,2)
hold on; box on;
plot(U_/Ud, imag(p_)/(2*pi), 'LineWidth', 2);
xlabel("U_{\infty}/U_D");
ylabel("p_I / 2\pi");

