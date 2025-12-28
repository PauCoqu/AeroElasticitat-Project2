function [modes,malla] = modal_analysis_wing(span, material,Nx,Ny,N_panels)
% Construim el model de viga i resolem el problema modal K*phi = w^2*M*phi.
% El resultat són les freqüències pròpies (f, omega) i les formes modals (Phi) de l'ala.

N_nodes = (Nx+1)*(Ny+1);
c = span.c_y;
y = span.y;
xi = linspace(0,1,Nx+1);  % Chordwise normalized coordinate (0 -> 1)


%% Malla ESTRUCTURAL 2D de viga:
   %*────*────*────*────*────*────* 
   %|    |    |    |    |    |    | 
   %*────*────*────*────*────*────*   -> y
   %|    |    |    |    |    |    | 
   %*────*────*────*────*────*────*  
   %y=0 y1   y2    y3   ...      y = b/2
   %* = node estructural

%Nodes
nodes = zeros(N_nodes,2); % 2 cops per -> [x,y]  
for i = 1:Nx+1 %i = direcció x
    for j = 1:Ny+1 %j direcció y
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

%%%%%%%%%%%%%%%%%%%%%%%%
figure; hold on; axis equal; grid on

% Elements
for e = 1:size(elem,1)
    xy = nodes(elem(e,:),:);
    plot([xy(:,1); xy(1,1)], [xy(:,2); xy(1,2)], 'k-');

    % número d’element (al centre)
    xc = mean(xy(:,1));
    yc = mean(xy(:,2));
    text(xc, yc, num2str(e), 'Color','r','FontSize',8,...
        'HorizontalAlignment','center');
end

% Nodes
plot(nodes(:,1), nodes(:,2), 'b.', 'MarkerSize', 12)

for k = 1:size(nodes,1)
    text(nodes(k,1), nodes(k,2), [' ' num2str(k)], ...
        'Color','b','FontSize',8);
end

xlabel('x'); ylabel('y');
title('Wing structural mesh (nodes + elements)');

%%%%%%%%%%%%%%%%%%%%%%%

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

%% Resolem el problema modal K_free * phi = w^2 * M_free * phi (Projecte 1)

[Phi, D] = eig(K_free, M_free);
omega = sqrt(diag(D));      % frequencies propies en rad/s
f     = omega/(2*pi);       % frequencies propies en Hz

[Phi_free, D] = eig(K_free, M_free);     % Phi_free: (#free DOF x #modes)
omega = sqrt(diag(D));
freq  = omega/(2*pi);

%Ordenem els modes segons frequencia creixent
[freq, idx] = sort(freq);
Phi_free = Phi_free(:,idx);

%Les 10 primeres freqüències (les més importants)
%fprintf('\n--- Modal analysis ---\n');
%disp('First 10 frequencies (Hz):');  
%disp(f_sorted(1:min(10, numel(f_sorted)))');

%% Guardem resultats
modes.Phi_free = Phi_free;   % size: (#free DOF) x nModes
modes.freq     = freq;
modes.M_free      = M_free;
modes.K_free      = K_free;
modes.freeDOFs = I_free;
modes.S        = S;
modes.Ix       = Ix;

malla.nodes = nodes;
malla.elements = elem;

end
