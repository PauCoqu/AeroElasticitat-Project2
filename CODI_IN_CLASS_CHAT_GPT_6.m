%% Control over a flat plate

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
lambda = 1;
c_root = c;
y0=0;

% Nelem
Nx = 10;
Ny = 20;

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

% Modes to keep
N_mod = 10;

% Gust params
sig_g = 1;
lambda_t = 762;
phi_gg = @(U_inf,w) sig_g^2*2*lambda_t./U_inf .* ...
    (1 + 8/3*(1.339*lambda_t*w./U_inf).^2) ./ ...
    (1 + (1.339*lambda_t*w./U_inf).^2).^(11/6);

%% Mesh
[nodes,elem,x_p,y_p,i_root,i_delta,i_e] = meshPlate(c,b,x0,y1,y2,xe,ye,Nx,Ny);

%% Structural matrices
[K,M,S,It,Ix,I_free] = structuralMatrices(nodes,elem,h,rho,E,nu,Nx,Ny,i_root);

% Number of DOFs
N_dof = size(K,1);

% Structural modes
[q_mod,f_mod] = structModes(K,M,I_free,N_mod);

%% Deflection vector kinematics
Id = zeros(N_dof,1);
for i = i_delta(:)'   % ensure row
    I_dof = 3*(i-1) + (1:3);
    Id(I_dof) = [x0-nodes(i,1); 0; 1];
end

%% Measurement
% Acceleration measurement vector
Ie = zeros(N_dof,1);
Ie(3*(i_e(1)-1)+1) = 1;

% Bending stress measurement
Se = zeros(N_dof,1);
Se(3*(i_e(2)-1)+2) = E*h/(2*nodes(i_e(2),2));

% Measurement tf function
%He = @(w) [-w^2*Ie', Se'];
He = @(w) [-w^2*Ie'; Se'];

%% Control law
Hd = [0.1, 0];
%Hd = [0, -100/(E*h)];
%Hd = [0.1, -100/(E*h)];

%% Loop through velocities and frequencies

%% Velocity vector
U_inf = 5:10:150; %desde 5 a 150 de 10 en 10

%% Frequency (Hz)
%f = [2, logspace(-5,0,5), 1.2:0.2:20];
f = [logspace(-5,0,5), 1:0.5:10];  %5 PUNTS DE 10^-5 FINS A 10^0
% log spacing for low f, linear spacing for high f

% Reduced matrices
K_red = q_mod'*K*q_mod;
M_red = q_mod'*M*q_mod;

% Initialize variables
Hg    = zeros(size(q_mod,2), length(f), length(U_inf));
Heg   = zeros(size(Hd,2),   length(f), length(U_inf));
Phi_gg = zeros(length(f), length(U_inf));
Phi_aa = zeros(length(f), length(U_inf));
Phi_ss = zeros(length(f), length(U_inf));

% Loop through velocities
for i = 1:length(U_inf)
    tic

    for j = 1:length(f)
        w = 2*pi*f(j);               % rad/s
        k = w*c/(2*U_inf(i));        % reduced freq

        % Aerodynamic matrices
        [Aeff,Ag,Ad] = aeroMatrices(rho_inf,a_inf,U_inf(i),k,c,x_p,y_p,S,Ix,It,Id,c_root,lambda,y0,b);

        % Reduced aero matrices
        A_red  = q_mod'*(Aeff + Ad*(Hd*He(w)))*q_mod;
        Ag_red = q_mod'*Ag;

        % System FRF
        Hg(:,j,i) = (-w^2*M_red + K_red - U_inf(i)^2*A_red) \ Ag_red;

        % Measurement FRF
        Heg(:,j,i) = He(w) * (q_mod*Hg(:,j,i));

        % Gust PSD
        Phi_gg(j,i) = phi_gg(U_inf(i),w);

        % Accel PSD
        Phi_aa(j,i) = abs(Heg(1,j,i))^2 * Phi_gg(j,i);

        % Stress PSD
        Phi_ss(j,i) = abs(Heg(2,j,i))^2 * Phi_gg(j,i);
    end

    fprintf('Velocity: %03i/%03i (U=%g m/s)\n', i, length(U_inf), U_inf(i));
end

% RMS values
g_ = squeeze(sqrt(trapz(f', Phi_gg, 1)));
a_ = squeeze(sqrt(trapz(f', Phi_aa, 1)));
s_ = squeeze(sqrt(trapz(f', Phi_ss, 1)));
disp([U_inf(:), a_(:), s_(:)])

%% PLOTS

% 1) RMS vs speed
figure; hold on; box on;
plot(U_inf, a_, 'LineWidth', 2);
xlabel('U_\infty [m/s]'); ylabel('a_{RMS} [m/s^2]');
title('Tip acceleration RMS vs U_\infty');

figure; hold on; box on;
plot(U_inf, s_, 'LineWidth', 2);
xlabel('U_\infty [m/s]'); ylabel('\sigma_{RMS} [Pa]');
title('Bending stress RMS vs U_\infty');

% 2) PSD curves for a few speeds
iu = round(linspace(1,length(U_inf), min(4,length(U_inf))));
figure; hold on; box on;
for kk = 1:length(iu)
    plot(f, Phi_aa(:,iu(kk)), 'LineWidth', 1.5);
end
set(gca,'XScale','log','YScale','log');
xlabel('f [Hz]'); ylabel('\Phi_{aa}');
title('Acceleration PSD');
legend(arrayfun(@(ii) sprintf('U=%.0f',U_inf(ii)), iu, 'UniformOutput', false), 'Location','best');

figure; hold on; box on;
for kk = 1:length(iu)
    plot(f, Phi_ss(:,iu(kk)), 'LineWidth', 1.5);
end
set(gca,'XScale','log','YScale','log');
xlabel('f [Hz]'); ylabel('\Phi_{\sigma\sigma}');
title('Stress PSD');
legend(arrayfun(@(ii) sprintf('U=%.0f',U_inf(ii)), iu, 'UniformOutput', false), 'Location','best');













%% Functions
% Mesh generation------------------------------------------------------------
function [nodes,elem,x_p,y_p,i_root,i_delta,i_e] = meshPlate(c,b,x0,y1,y2,xe,ye,Nx,Ny)
% (based on class mesh logic)

Nx2 = ceil(Nx*(c-x0)/c);
Nx1 = Nx - Nx2;

Ny3 = ceil(Ny*(b-y2)/b);
Ny2 = ceil(Ny*(y2-y1)/b);
Ny1 = Ny - Ny2 - Ny3;

% x coordinates (avoid duplicate points at joins)
x1 = linspace(0,  x0, Nx1+1);
x2 = linspace(x0, c,  Nx2+1); x2 = x2(2:end);
x  = [x1, x2];

% y coordinates (avoid duplicate points at joins)
y1v = linspace(0,  y1, Ny1+1);
y2v = linspace(y1, y2, Ny2+1); y2v = y2v(2:end);
y3v = linspace(y2, b,  Ny3+1); y3v = y3v(2:end);
y   = [y1v, y2v, y3v];

% Nodes
nodes = zeros((Nx+1)*(Ny+1),2);
for ix = 1:(Nx+1)
    for iy = 1:(Ny+1)
        k = (Ny+1)*(ix-1) + iy;
        nodes(k,1) = x(ix);
        nodes(k,2) = y(iy);
    end
end

% Elements and control-surface marking (per element)
elem = zeros(Nx*Ny,4);
isCtrl = false(Nx*Ny,1);

for ix = 1:Nx
    for iy = 1:Ny
        e = Ny*(ix-1) + iy;

        elem(e,1) = (Ny+1)*(ix-1) + iy;
        elem(e,2) = (Ny+1)* ix    + iy;
        elem(e,3) = (Ny+1)* ix    + iy + 1;
        elem(e,4) = (Ny+1)*(ix-1) + iy + 1;

        % Control surface: x in (x0..c) -> ix > Nx1, and y in (y1..y2)
        if (ix > Nx1) && (iy > Ny1) && (iy <= Ny1 + Ny2)
            isCtrl(e) = true;
        end
    end
end

% Control surface nodes
i_delta = unique(elem(isCtrl,:));

%% AERODYNAMIC MESH (panel points)
x_p = zeros(Nx*Ny,3);
y_p = zeros(Nx*Ny,3);

for e = 1:Nx*Ny
    % quarter-chord points on each side + collocation
    x_p(e,1) = nodes(elem(e,1),1) + (nodes(elem(e,2),1) - nodes(elem(e,1),1))/4;
    y_p(e,1) = nodes(elem(e,1),2);

    x_p(e,2) = nodes(elem(e,4),1) + (nodes(elem(e,3),1) - nodes(elem(e,4),1))/4;
    y_p(e,2) = nodes(elem(e,4),2);

    x_p(e,3) = nodes(elem(e,1),1) + 3*(nodes(elem(e,2),1) - nodes(elem(e,1),1))/4;
    y_p(e,3) = (nodes(elem(e,1),2) + nodes(elem(e,4),2))/2;
end

% Root nodes (x=0 line)
i_root = 1:(Ny+1):size(nodes,1);

% Closest node(s) to sensors
i_e = zeros(length(xe),1);
for k = 1:length(xe)
    [~, i_e(k)] = min(abs(nodes(:,1) - xe(k)) + abs(nodes(:,2) - ye(k)));
end
end

% Structural matrices computation ------------------------------------------
function [K,M,S,It,Ix,I_free] = structuralMatrices(nodes,elem,h,rho,E,nu,Nx,Ny,i_root)

N_dof = 3*size(nodes,1);
M = sparse(N_dof,N_dof);
K = sparse(N_dof,N_dof);
S = sparse(N_dof,Nx*Ny);
Ix = sparse(Nx*Ny,N_dof);
It = sparse(Nx*Ny,N_dof);

for e = 1:Nx*Ny
    % Element size
    a_e = (nodes(elem(e,2),1) - nodes(elem(e,1),1))/2;
    b_e = (nodes(elem(e,4),2) - nodes(elem(e,1),2))/2;

    % Element matrices (must exist in your path)
    M_e = plateMass(a_e,b_e,h,rho);
    K_e = plateStiffness(a_e,b_e,h,E,nu);
    S_e = plateForce(a_e,b_e);

    % Interpolation (force these as 1x12 rows)
    Ix_e = reshape([
        -9/(32*a_e), -9*b_e/(64*a_e),  5/32, ...
         9/(32*a_e),  9*b_e/(64*a_e), -5/32, ...
         9/(32*a_e), -9*b_e/(64*a_e),  5/32, ...
        -9/(32*a_e),  9*b_e/(64*a_e), -5/32 ...
    ], 1, 12);

    It_e = reshape([
         5/64,  5*b_e/128, -3*a_e/64, ...
        27/64, 27*b_e/128,  9*a_e/64, ...
        27/64, -27*b_e/128, 9*a_e/64, ...
         5/64, -5*b_e/128, -3*a_e/64 ...
    ], 1, 12);

    % Global DOF indices (12x1)
    I_dof = [
        3*(elem(e,1)-1) + (1:3);
        3*(elem(e,2)-1) + (1:3);
        3*(elem(e,3)-1) + (1:3);
        3*(elem(e,4)-1) + (1:3);
    ];
    I_dof = I_dof(:);

    % Assembly
    M(I_dof,I_dof) = M(I_dof,I_dof) + M_e;
    K(I_dof,I_dof) = K(I_dof,I_dof) + K_e;
    S(I_dof,e)     = S(I_dof,e)     + S_e;

    Ix(e,I_dof) = Ix(e,I_dof) + Ix_e;
    It(e,I_dof) = It(e,I_dof) + It_e;
end

% Clamp root
I_fix = [3*(i_root-1)+1, 3*(i_root-1)+2, 3*(i_root-1)+3];
I_fix = I_fix(:);
I_free = setdiff(1:N_dof, I_fix);
end

% Structural modes ---------------------------------------------------------
function [q_mod,f_mod] = structModes(K,M,I_free,N_mod)

N_dof = size(K,1);
q_mod = zeros(N_dof, N_mod);

% Eigs on free DOFs
[q_free, w2] = eigs(K(I_free,I_free), M(I_free,I_free), N_mod, 'smallestabs');

q_mod(I_free,:) = q_free;
f_mod = sqrt(diag(w2))/(2*pi);
end

% Aerodynamic matrices computation ----------------------------------------
function [Aeff,Ag,Ad] = aeroMatrices(rho_inf,a_inf,U_inf,k,c,x_p,y_p,S,Ix,It,Id,c_root,lambda,y0,b)

M_inf = U_inf/a_inf;

% AIC matrix (must exist, else fallback stub below)
AIC = buildAIC(x_p,y_p,c_root,lambda,y0,b,M_inf,k);
%buildAIC(x_p,y_p,k,M_inf,c/2);

% Gust vector delay
Ig = exp(1i*k*2/c*x_p(:,3));

Aeff = rho_inf/2 * S * (AIC \ (1i*k*2/c*It + Ix));
Ag   = -rho_inf/2 * S * (AIC \ Ig);
Ad   = Aeff * Id;
end



