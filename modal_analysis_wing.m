function [modes] = modal_analysis_wing(span, material)
% Construim el model de viga i resolem el problema modal K*phi = w^2*M*phi.
% El resultat són les freqüències pròpies (f, omega) i les formes modals (Phi) de l'ala.

y  = span.y;    % (semiala)
EI = span.EI;   
GJ = span.GJ;   
m  = span.m;
E  = material.YoungModulus;     
nu = material.Poisson;          
G  = E/(2*(1+nu));              

%% Model 1D de viga en y:
   %────*────*────*────*────*────* -> y
   %y=0 y1  y2   y3   ...       y = b/2
   %* = node estructural

Ne   = numel(y) - 1; % nombre d'elements
Nn   = numel(y); % nombre de nodes
nd   = 3; % DOFs per node: [eta (desplaç.), zeta (rot.), theta (torsió)]
ndof = nd * Nn; % nombre total de graus de llibertat

% Matrius globals de rigidesa i massa
K = sparse(ndof, ndof);
M = sparse(ndof, ndof);

%% Ensamblatje
for e = 1:Ne   % Índex dels dos nodes de l'element e
    i1 = e;       % node inicial de l'element
    i2 = e + 1;   % node final de l'element
    L = y(i2) - y(i1); % Longitud de l'element

    % Propietats mitjanes de l'element (promig)
    EIe = 0.5*(EI(i1) + EI(i2));   
    GJe = 0.5*(GJ(i1) + GJ(i2));   
    me  = 0.5*(m(i1)  + m(i2));    % massa per unitat de longitud
    Ie = EIe / E;   % moment d'inèrcia a flexió
    Je = GJe / G;   % constant de torsió
    I_sc = 0.5*(span.I_sc(i1) + span.I_sc(i2));
    x_cm = 0.5*(span.x_cm(i1) + span.x_cm(i2));
    x_sc = 0.5*(span.x_sc(i1) + span.x_sc(i2));

    % Matrius elementals de massa i rigidesa
    Me = beamMass(L, me, I_sc, x_cm, x_sc);
    Ke = beamStiffness(L, E, Ie, G, Je);

    % Calculem els DOFs associats a aquest element (3 DOFs per node)
    %Formula que funciona
    %(i1-1)*nd + (1:nd)   = (1-1)*3 + [1 2 3] = [1 2 3]
    %(i2-1)*nd + (1:nd)   = (2-1)*3 + [1 2 3] = [4 5 6]
    dofs = [ (i1-1)*nd + (1:nd), (i2-1)*nd + (1:nd) ];
      
    % Sumem la contribució de l'element a la matriu global
    M(dofs, dofs) = M(dofs, dofs) + Me;
    K(dofs, dofs) = K(dofs, dofs) + Ke;
end

%% Boundary conditions: Arrel encastada
% Fixem els tres DOFs del primer node (arrel de l'ala):  eta(1) = 0, zeta(1) = 0, theta(1) = 0
fixed = [1 2 3]; % DOFs bloquejats

%del 1:nof extreiem els fixed values
free  = setdiff(1:ndof, fixed);   % setdiff treu elements d'un conjunt.
Kff = K(free, free);
Mff = M(free, free);

%Estudiem tots els DOF exepte els 3 primers.

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
