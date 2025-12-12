% Control over a flat plate

clear
close all
clc

%% Input data
% Plate geometrical params
c = 0.75;
b = 2;
h = 0.01;
x0 = 0.7*c;
y1 = 0.5*b;
y2 = 0.9*b;

% Nelem
Nx = 10;
Ny = 20;
N_nod = 6;

% Sensor location
xe = [0.5*c, 0.5*c];
ye = [0.75*b, b/Ny];

% Material props
E = 70e9;
nu = 0.3;
rho = 2300;

% Flow properties
rho_inf = 1.25;
a_inf = 343;


% Gust params
sig_g = 1;
lambda_t = 762;
phi_gg = @(U_inf,w) sig_g^2*2*lambda_t./U_inf*(1+8/3*(1.339*lambda_t*w./U_inf).^2)./(1+(1.339*lambda_t*w./U_inf).^2).^(11/6);

%% Mesh
[nodes,elem,x_p,y_p,i_root,i_delta,i_e] = meshPlate(c,b,x0,y1,y2,xe,ye,Nx,Ny);

%% Structural matrices
[K,M,S,It,Ix,I_free] = structuralMatrices(nodes,elem,h,rho,E,nu,Nx,Ny);

% Number of DOFs
N_dof = size(K,1);

% Structural modes
[q_mod,f_mod] = structModes(K,M,I_free,N_mod);

% Deflection vector kinematics
Id = zeros(N_dof,1);
for i = i_delta'
    I_dof = 3*(i-1) + [1,2,3];
    Id(I_dof) = [
        x0-nodes(i,1);
        0;
        1;
        ];
end

%% Measurement
% Acceleration measurement vector
Ie = zeros(N_dof,1);
Ie(3*(i_e(1)-1)+1) = 1;

% Bending stress measurement
Se = zeros(N_dof,1);
Se(3*(i_e(2)-1)+2) = E*h/(2*nodes(i_e(2),2));

% Measurement tf function
He = @(w) [-w^2*Ie',Se'];

%% Control law
% Control law matrix
Hd = [0.1,0];
%Hd = [0, -100/(E*h)];
%Hd = [0.1, -100/(E*h)];

%% Loop through velocities and frequencies

% Velocity vector
U_inf = 20:5:180;

%Frequency
f = [0, logspace(-5,0,50),1.2:0.2:20];

% Reduced matrices
K_red = q_mod'*K*q_mod;
M_red = q_mod'*M*q_mod;

% Initializa variables
Hg = zeros(size(q_mod,2),length(f),length(U_inf));
Heg = zeros(size(Hd,2),length(f),length(U_inf));
Phi_gg = zeros(length(f),length(U_inf));
Phi_aa = zeros(length(f),length(U_inf));
Phi_ss = zeros(length(f),length(U_inf));

% Loop through velocities
for i = 1:length(U_inf)
    %Init timer
    tic

    % Loop through freqs
    for j = 1:length(f)
        %rad/S
        w = 2*pi*f(j);

        % Reduced freq
        k = w*c/(2*U_inf(i));

        %Aerodynamic mat
        [Aeff,Ag,Ad] = aeroMatrices(rho_inf,a_inf,U_inf(i),k,c,x_p,y_p,S,Ix,It,Id);

        %Reduced aero mat
        A_red = q_mod'*(Aeff + Ad*(Hd*He(w)))*q_mod;
        Ag_red = q_mod'*Ag;

        %System FRF
        Hg(:,j,i) = (-w^2*M_red + K_red - U_inf(i)^2*A_red)\Ag_red;

        %Measurement FRF 
        Heg(:,j,i) = He(w)*(q_mod*Hg(:,j,i));

        % Gust PSD
        Phi_gg(j,i) = phi_gg(U_inf(i),w);

        %Accel PSD
        Phi_aa(j,i) = abs(Heg(1,j,i))^2*Phi_gg(j,i); 

        %Stress PSD
        Phi_ss(j,i) = abs(Heg(2,j,i))^2*Phi_gg(j,i); 
    end

    % Show progress
    fprintf('Velocity: %03i/%03i\n',i,length(U_inf));
end


%Get RMS values
g_ = squeeze(sqrt(trapz(f',Phi_gg,1)));
a_ = squeeze(sqrt(trapz(f',Phi_aa,1)));
s_ = squeeze(sqrt(trapz(f',Phi_ss,1)));
%% Functions
% Mesh generation------------------------------------------------------------
function [nodes,elem,x_p,y_p,i_root,i_delta,i_e] = meshPlate(c,b,x0,y1,y2,xe,ye,Nx,Ny)
% COPIED FROM CLASS 3
Nx2 = ceil(Nx*(c-x0)/c); % Changed xh to x0!
Nx1 = Nx - Nx2;

Ny3 = ceil(Ny*(b-y2)/b);
Ny2 = ceil(Ny*(y2-y1)/b);
Ny1 = Ny - Ny2 - Ny3;

x = [0:x0/Nx1:x0,x0+(c-x0)/Nx2:(c-x0)/Nx2:c]; %Crec que aixo esta mal

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

% Get the root nodes
i_root = 1:Ny+1:size(nodes,1);

% Closest node to sensor
i_e = zeros(length(xe),1);
for i = 1:length(xe)
    [~,i_e(i)] = min(abs(nodes(:,1) - xe(i)) + abs(nodes(:,2) - ye(i)));
end

end


% Structural matrices computation ------------------------------------------
function [K,M,S,It,Ix,I_free] = structuralMatrices(nodes,elem,h,rho,E,nu,Nx,Ny,i_root)
    % COPIED FROM CLASS 4/5
% Assembly of the structural matrices
%Assembly of structural matrices
N_dof = 3*size(nodes,1);
M = sparse(N_dof,N_dof);
K = sparse(N_dof,N_dof);
S = sparse(N_dof,Nx*Ny);
Ix = sparse(Nx*Ny,N_dof);
It = sparse(Nx*Ny,N_dof);

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
        ];

    % Collocation point interpolation
    It_e = [
        5/64, 5*b_e/128, -3*a_e/64,
        27/64, 27*b_e/128, 9*a_e/64,
        27/64, -27*b_e/128, 9*a_e/64,
        5/64, -5*b_e/128, -3*a_e/64,
        ];


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

I_fix = 3*(repmat(i_root,3,1)-1)+[1,2,3];
I_free = setdiff(1:N_dof,I_fix);
end

% Structural modes -------------------------------------------------------
function [q_mod,f_mod] = structModes(K,M,I_free,N_mod)

%Number of dofs
N_dof = size(K,1);

%Vibration modes
q_mod = zeros(N_dof,N_mod);
[q_mod(I_free,:),w2] = eigs(K(I_free,I_freee),M(I_free,I_free),N_mod,'smallestabs');

% Natural freqs
f_mod = sqrt(diag(w2))/(2*pi);


end



% Aerodynamic matrices computation ---------------------------------------
function [Aeff,Ag,Ad] = aeroMatrices(rho_inf,a_inf,U_inf,k,c,x_p,y_p,S,Ix,It,Id)

% Mach number
M_inf = U_inf/a_inf;

% AIC matrix
AIC = dlmAIC_vec(x_p,y_p,k,M_inf,c/2);

%Gust vector delay
Ig = exp(1i*k*2/c*x_p(:,3));

%Aeroelastic coupling matrices
Aeff = rho_inf/2*S*(AIC\(1i*k*2/c*It + Ix));

% Gust loads vector
Ag = -rho_inf/2*S*(AIC\Ig);

% Deflection load vector
Ad = Aeff*Id;
end
