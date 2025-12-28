function [modes] = modal_analysis_wing(span, material,Nx,Ny)
% Construim el model de viga i resolem el problema modal K*phi = w^2*M*phi.
% El resultat són les freqüències pròpies (f, omega) i les formes modals (Phi) de l'ala.

N_panels = Nx*Ny;
N_nodes = (Nx+1)*(Ny+1);
x_p = zeros(N_panels,3);
y_p = zeros(N_panels,3);
c = span.c_y;
y = span.y;
b = span.b;
h = span.c_y*0.05; %en cas de NACA 0010

xi = linspace(0,1,Nx+1);  % Chordwise normalized coordinate (0 -> 1)


%% Model ESTRUCTURAL 2D de viga:
   %*────*────*────*────*────*────* 
   %|    |    |    |    |    |    | 
   %*────*────*────*────*────*────*   -> y
   %|    |    |    |    |    |    | 
   %*────*────*────*────*────*────*  
   %y=0 y1   y2    y3   ...      y = b/2
   %* = node estructural

%Nodes
nodes = zeros(N_nodes,2); % 2 cops per -> [x,y]  
for i = 1:Nx+1
    for j = 1:Ny+1
        k = (Ny+1)*(i-1) + j;
        nodes(k,1) = xi(i) * c(j);  %EXPLICAR MOTIU NORMALITZACIÓ
        nodes(k,2) = y(j);
    end
end

%Elements
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

%% MALLA AERODINAMICA

for e = 1:N_panels
    n1 = elem(e,1); 
    n2 = elem(e,2); 
    n3 = elem(e,3); 
    n4 = elem(e,4);

    % quarter-chord endpoints (points 1 and 2)
    x_p(e,1) = nodes(n1,1) + (nodes(n2,1)-nodes(n1,1))/4;
    y_p(e,1) = nodes(n1,2);

    x_p(e,2) = nodes(n4,1) + (nodes(n3,1)-nodes(n4,1))/4;
    y_p(e,2) = nodes(n4,2);

    % collocation (3/4 chord, mid-span of element)
    x_p(e,3) = nodes(n1,1) + 3*(nodes(n2,1)-nodes(n1,1))/4;
    y_p(e,3) = 0.5*(nodes(n1,2)+nodes(n4,2));
end

%% Ensamblatje matrius estructurals

N_dof = 3*N_nodes;
M = zeros(N_dof,N_dof);
K = zeros(N_dof,N_dof);
S = zeros(N_dof,N_panels);
Ix = zeros(N_panels,N_dof);

for e = 1:N_panels    
    % element half-sizes
    b_e = 0.5*(nodes(elem(e,4),2) - nodes(elem(e,1),2)); %constant
    
    % a_e: és variable quan tenim tapper. (0.25 perque dividim per 2 dos cops)
    
    a_e = 0.25*((nodes(elem(e,2),1) - nodes(elem(e,1),1)) + (nodes(elem(e,3),1) - nodes(elem(e,4),1))) ;
   
    %mid chord
    y_mid = 0.5*(nodes(elem(e,1),2) + nodes(elem(e,4),2));
    c_mid = span.c_root*(1-(1-span.lambda)*(y_mid/span.b));
    h_e = 0.10 * c_mid;

    M_e = plateMass(a_e,b_e,h_e,material.Density);
    K_e = plateStiffness(a_e,b_e,h_e,material.YoungModulus,material.Poisson);
    S_e = plateForce(a_e,b_e);

    %Coeficients del professor
    Ix_e = [-9/(32*a_e), -(9*b_e)/(64*a_e), 5/32, 9/(32*a_e), (9*b_e)/(64*a_e), -5/32, 9/(32*a_e), -(9*b_e)/(64*a_e), 5/32, -9/(32*a_e), (9*b_e)/(64*a_e), -5/32 ];

    % global DOF indices for the 4 nodes
    I_dof = [3*(elem(e,1)-1) + (1:3)',
             3*(elem(e,2)-1) + (1:3)', 
             3*(elem(e,3)-1) + (1:3)',
             3*(elem(e,4)-1) + (1:3)' ];

    % assembly
    M(I_dof,I_dof) = M(I_dof,I_dof) + M_e;
    K(I_dof,I_dof) = K(I_dof,I_dof) + K_e;
    S(I_dof,e)     = S(I_dof,e) + S_e;
    Ix(e,I_dof)    = Ix(e,I_dof) + Ix_e;
end

%% Boundary conditions: Arrel encastada
% Fixem els tres DOFs de tots els nodes en y=0 (arrel de l'ala):  [eta (desplaç.)=0, zeta (rot.)=0, theta (torsió)=0]
I_fix = [
   1:3*(Ny+1):N_dof,...
   2:3*(Ny+1):N_dof,...
   3:3*(Ny+1):N_dof
   ];
I_free = setdiff(1:N_dof,I_fix); % setdiff treu elements d'un conjunt.

%Estudiem tots els DOF exepte aquests

%Reduced matrices
M_free = M(I_free,I_free);
K_free = K(I_free,I_free);
S_free = S(I_free,:);

nd   = 3; % DOFs per node: [eta (desplaç.), zeta (rot.), theta (torsió)]
ndof = nd * N_nodes; % nombre total de graus de llibertat

%%%%%%%%%%%%%%%%%%%%%%%%

% Aeroelastic coupling
rho_inf = 1.25;
A0 = rho_inf/2*S*(AIC\Ix);

% Comput divergence speed
[X_div,U_inf] = eigs(sparse(K(I_free,I_free)),sparse(A0(I_free,I_free)),10,"sm");

%%%%%%%%%%%%%%%%%%%%%%%%%


%% Resolem el problema modal Kff*phi = w^2*Mff*phi

[Phi, D] = eig(full(Kff), full(Mff));
omega = sqrt(diag(D));      % frequencies propies en rad/s
f     = omega/(2*pi);       % frequencies propies en Hz

% Ordenem els modes segons frequencia creixent
[f_sorted, idx] = sort(f);
Phi = Phi(:, idx);

%Les 10 primeres freqüències (les més importants)
fprintf('\n--- Modal analysis ---\n');
disp('First 10 frequencies (Hz):');  
disp(f_sorted(1:min(10, numel(f_sorted)))');

%% Guardem resultats
modes.f        = f_sorted;    % frequencies propies [Hz]
modes.omega    = omega(idx);  % frequencies propies [rad/s]
modes.Phi      = Phi;         % modes propis (només DOFs lliures)
modes.freeDOFs = free;        % index dels DOFs lliures
modes.Mff      = Mff;         % matriu de massa reduïda
modes.Kff      = Kff;         % matriu de rigidesa reduïda

end
