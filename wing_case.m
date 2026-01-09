function metrics = wing_case(density,youngModulus,poissonRatio, b, lambda, Nx, Ny,U_inf,y0)

%% STEP 1: Section properties (del pefil) ('beam_properties')

% Section properties
geometry = "naca"; %Perfil Naca (ENUNCIAT)
param = [0.85,0.085];  %param(1)= chord (ENUNCIAT); param(2) = thickness (0.1*chord per al NACA 0010)
resolution = 2; %Canvia el mallat del perfil (Figure 1)  
nodesFile = ""; %EN BLANC (es per crear una malla des d'un fitxer)
elemsFile = ""; %EN BLANC (es per crear una malla des d'un fitxer)

% Obtain mesh data (true ploteja el grafic, false no)
[nodes,material,elems,elemat] = getMeshData(geometry,param,resolution,density,youngModulus,poissonRatio,nodesFile,elemsFile,false);

c_root = param(1); %corda root
material.YoungModulus = youngModulus;
material.Poisson      = poissonRatio;
material.Density      = density;

% Compute section properties (%true ploteja el grafic, false no)
section = getSectionProperties(nodes,material,elems,elemat,false);

% The warping displacement:
scale =0; % de 0 a 1, deflexió del perfil

%% STEP 2: Spanwise distributions 
    
ctip = lambda * c_root; %corda punta
Surf    = b * (c_root + ctip); % Area de les DOS ales
AR   = (2*b)^2 / Surf; % Aspect ratio 

% Corda a cada punt: c(y) = c_root*[1 - (1-lambda)*2y/b]
y = linspace(y0, b+y0, Ny+1)'; % coordenades dels nodes tenint en compte y0=0.34
c_y = c_root * (1 - (1 - lambda) * (y-y0)/(b)); %corda cada perfil

%Calcul de les propietats al llarg de l'envergadura [span]
[span] = spanwise_distributions(section, material, c_y, c_root, y, b, Surf, AR, lambda);

%% STEP 3 - Kirchoff_love_plate 

% AQUEST PAS NOMÉS CAL EXECUTAR-LO UNA VEGADA!!!
%Cal correr kirchoff_love_plate
%En la funció (kirchoff_love_plate) es genera les funcions plateMass.m, plateStiffness.m y plateForce.m.
% Després al modal analysis extreiem les matrius: M_e (masa) ,K_e (rigidesa),B_e (cargues aerodinamiques)

%% STEP 4 - Malla ESTRUCTURAL i Modal Analysis 

N_panels = Nx*Ny;
[modes,malla,N_dof] = modal_analysis_wing(span, material,Nx,Ny,N_panels); %Calculem les matrius
%f1 = modes.f(1);   % primera freq pròpia

%% STEP 5 - Malla AERODINAMICA

x_p = zeros(N_panels,3);
y_p = zeros(N_panels,3);

nodes = malla.nodes;
elem = malla.elements;

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

%% STEP 5.1 - Aerodynamic influence coefficients

% Fluid data
rho_inf = 1.25;
%a_inf = 343;
%eta = 0.1*(c/2);  %Amplitud del moviment vertical (CAS ESTATIC)
%freq = 1; % Hz (CAS ESTATIC)
%alpha = alpha_deg*pi/180;  %rad

% Reduced frequency -> Frequencia d'oscilació (k>0 -> flutter i altres, k=0 -> CAS ESTÀTIC)
k = 0; %2*pi*freq*c/(2*U_inf)
M_inf = 0; % Mach = U_inf/a_inf; %0; %CAS INCOMPRESSIBLE;

% Wref = alpha*ones(N_panels,1); PROJECTE 1
%Wref = -1i*k*eta/(c/2)*ones(N_panels,1); 

% AIC matrix coefficients
AIC = zeros(N_panels,N_panels);
for i = 1:N_panels
    % Collocation point coordinates of panel "i"
    x_i = x_p(i,3);
    y_i = y_p(i,3);
    % Loop through all doublet segments
    for j = 1:N_panels
        y_mid_j = y_p(j,3);
        c_local = c_root * (1 - (1 - lambda) * (y_mid_j - y0)/b);
        % We add the induced velocity contribution to the AIC
        AIC(i,j) = AIC(i,j) + w_doublet(x_i,y_i,x_p(j,:),y_p(j,:),c_local/2,M_inf,k,1);
        AIC(i,j) = AIC(i,j) + w_doublet(x_i,y_i,x_p(j,[2,1,3]),-y_p(j,[2,1,3]),c_local/2,M_inf,k,1);
    end
end

%% 1. DIVERGENCE SPEEED
% It can be done with the DOFs or with an approximation onto some structural modes (which need modal_analysis_wing.m)
% Ho farem amb tots els DOFs primer

%Redefinim totes les variables guardades al modal_analysis_wing
M       = modes.M;
K       = modes.K;
S       = modes.S;
Ix      = modes.Ix;
It      = modes.It;
K_free  = modes.K_free;
M_free  = modes.M_free;
I_free  = modes.freeDOFs;
S_free  = modes.S_free;
Phi_free= modes.Phi_free;


A0 = 0.5*rho_inf * (S*(AIC\Ix)); %This condition is the resultant of Wref = -Uinf*Ix*q when k=0????
A0_free = A0(I_free,I_free);

% Compute divergence speed(s): Kff x = U^2 A0ff x ????
%it extract the eigenvectors and eigenvalues. We only want the eigenvalues U^2
[~,U2] = eigs(sparse(K_free), sparse(A0_free), 10, "smallestabs");
U = sort(sqrt(real(diag(U2)))); 
Ud = U(1);



%% 2. RISK OF FLUTTER (ANEM PER AQUÍ ARA)
% for i=1:length(U_inf)
%     Wref = -U_inf(i) * (Ix .* q_free);
% end
%Wref = 1; 

%--------------------------
%% p method

%k = 0; % quasi-steady  %JA DEFINIT MÉS AMUNT
%M_inf = 0; % Incompressibility %JA DEFINIT MÉS AMUNT

%Reuse modal results from modal_analysis_wing
q_mod = zeros(N_dof, size(modes.Phi_free,2));
q_mod(I_free,:) = Phi_free;

% Select modes for p-method as we do a Model reduction
i_modes = [1,2,3,5];
i_modes = [1,2];  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BORRAR
N_modes = length(i_modes);
Phi = q_mod(:, i_modes);

% Reduced matrices
M_red  = Phi.'*M*Phi;
K_red  = Phi.'*K*Phi;
S_red  = Phi.'*S;
Ix_red = Ix*Phi;
It_red = It*Phi;

A0_red = rho_inf/2*S_red*(AIC\Ix_red);
A1_red = -rho_inf/2*S_red*(AIC\It_red);

% Initialize velocity vector
U_ = 5:5:1.4*Ud;

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

    %Check stability
    if i > 1 && max(real(p_(:,i))) > tol && isempty(U_min) % Check if the min vel is encounterd. If yes, the flutter speed is betweem this iteration and the past one
    U_min = U_(i-1);
    U_max = U_(i);
    break;
    end
end

% PLOTS

figure
subplot(2,1,1)
hold on; box on;
plot(U_./Ud,real(p_).*c_root./(2*U_));
xlabel("U_{\infty}/U_D");
ylabel("p_Rc/2U");

subplot(2,1,2)
hold on; box on;
plot(U_./Ud,imag(p_)./(2*pi));
xlabel("U_{\infty}/U_D");
ylabel("p_I/2\pi");

%--------------------------

% %% k method
% % We need to recompute the matrices at each iteration
% 
% % Aerodynamic
% M_inf = 0;
% AIC = @(k_) dlmAIC_vec(x_p,y_p,k_,M_inf,c/2);
% 
% % Initialized reduced freq
% dinvk = 0.1;
% invk_ = dinvk:dinvk:20;
% 
% % Zero tolerance
% tol = 1e-6;
% 
% % Initialize variables
% m_ = nan(N,N,length(invk_));
% l_ = nan(N,length(invk_));
% Umin = [];
% Umax = [];
% 
% % Loop through reduced frequencues
% for i = 1:length(invk_)
% 
%     % Initilize timer
%     tic
% 
% 
% 
%     % Effective matrices
%     k = 1/invk-(i);
%     Beff = M_red + rho_inf/2*S_red*(AIC(k)\(1i*c/(2*k)*It_red + c^2/(4*k^2)*Ix_red));
% 
%     %Eigenvalues
%     [X, L] = eig(K_red,Beff);
%     l = 1./diag(L);
% 
% 
%     % Sort modes
%     if i == 1
%         l_(:,i) = l;
%         m_(:,:,i) = X;
%     else
%         m2sort = 1:N;
%         for j = 1:length(l)
%             [~,jast] = min(abs(real(l_(j,i-1))-real(l(m2sort))) + abs(imag(l_(j,i-1))-imag(l(m2sort))));
%             l_(j,i) = l(m2sort(jast));
%             m_(:,j,i) = X(:,m2sort(jast));
%             m2sort = setdiff(m2sort,m2sort(jast));
%         end
%     end
% 
%     %Check stability
%     if max(imag(l_(:,i))) > tol && isempty(Umin)
% 
%         [~,jast] = max(imag(l_(:,i)));
%         Umin = sqrt(1./real(l_(jmax,i-1)))*c(2*invk_(i-1));
%         Umax = sqrt(1./real(l_(jmax,i)))*c(2*invk_(i-1));
%     end
% 
%     % Print iteration time
%     %ja tenim el toc, ha jo fare aixo
% end
% 
% 
% % Recover parameters
% g_ = imag(l_)./real(l_);
% w_ = sqrt(1./real(l_));
% U_ = w_*c/2.*invk_;
% 
% % Plots
% 
%     figure(2)
%     subplot(2,1,1)
%     hold on; box on;
%     plot(U_'/Ud,g_');
%     xlabel("U/U_D");
%     ylabel("g");
% 
%     subplot(2,1,2)
%     hold on; box on;
%     plot(U_'/Ud,w_'/(2*pi));
%     xlabel("U/U_D");
%     ylabel("w/2pi");

%% k method (DLM) — based on your class notes

% --------------------------
% Aerodynamic AIC as function of k
M_inf = 0;
c_ref = 0.5*(c_root + ctip);

AICfun = @(k_) AIC_builder(k_, x_p, y_p, N_panels, c_ref, lambda, y0, b, M_inf);


% Sweep in invk = 1/k
dinvk = 0.1;
invk_ = dinvk:dinvk:20;        % invk = 1/k
Nk    = length(invk_);

tol = 1e-6;

% Storage
N = size(K_red,1);             % number of retained modes
m_ = nan(N,N,Nk);              % eigenvectors (modal coordinates)
l_ = nan(N,Nk);                % store lambda' = 1/lambda (complex)

Umin = [];
Umax = [];
jmax = [];

% Loop through reduced frequencies
for i = 1:Nk

    k = 1/invk_(i);

    % Build AIC(k)
    AICk = AICfun(k);

    % Beff(k)
    % Beff = M_red + (rho/2) * S_red * AIC^{-1} * ( i*c/(2k)*It_red + c^2/(4k^2)*Ix_red )
    AUX  = (1i*c_ref/(2*k))*It_red + (c_ref^2/(4*k^2))*Ix_red;  % per no fallar parentesis
    Beff = M_red + (rho_inf/2) * (S_red * (AICk \ AUX));

    % Generalized eigenproblem: K_red * x = lambda * Beff * x
    [X, L] = eig(K_red, Beff);
    lam = diag(L);

    % Store lambda' = 1/lambda
    l = 1./lam;

    % Mode tracking
    if i == 1
        l_(:,i)   = l;
        m_(:,:,i) = X;
    else
        m2sort = 1:N;
        for j = 1:N
            % nearest in complex plane
            [~, jast] = min( abs(real(l_(j,i-1)) - real(l(m2sort))) + ...
                             abs(imag(l_(j,i-1)) - imag(l(m2sort))) );
            l_(j,i)   = l(m2sort(jast));
            m_(:,j,i) = X(:,m2sort(jast));
            m2sort    = setdiff(m2sort, m2sort(jast));
        end
    end

    % Check flutter bracket using g crossing (more correct than imag(l)>0)
    if i > 1 && isempty(Umin)
        g_prev = imag(l_(:,i-1)) ./ real(l_(:,i-1));
        g_curr = imag(l_(:,i))   ./ real(l_(:,i));

        % find a mode that crosses 0
        cross = find(g_prev .* g_curr < 0, 1, "first");
        if ~isempty(cross)
            jmax = cross;

            % Recover omega and U at the two points (from lambda')
            w_prev = sqrt(1./real(l_(jmax,i-1)));
            w_curr = sqrt(1./real(l_(jmax,i)));

            U_prev = w_prev * c_ref/2 * invk_(i-1);
            U_curr = w_curr * c_ref/2 * invk_(i);

            Umin = U_prev;
            Umax = U_curr;
        end
    end

%     %Check stability
%     if max(imag(l_(:,i))) > tol && isempty(Umin)
% 
%         [~,jast] = max(imag(l_(:,i)));
%         Umin = sqrt(1./real(l_(jmax,i-1)))*c(2*invk_(i-1));
%         Umax = sqrt(1./real(l_(jmax,i)))*c(2*invk_(i-1));
%     end
% 
end

% Recover parameters
g_ = imag(l_)./real(l_);
w_ = sqrt(1./real(l_));          % omega
U_ = (w_*c_ref/2).*invk_;      % U = omega*c/2 * invk


% Plots
figure
subplot(2,1,1)
hold on; box on;
plot(U_.'/Ud, g_.', 'LineWidth', 1.2);
yline(0,'--');
xlabel("U/U_D");
ylabel("g");

subplot(2,1,2)
hold on; box on;
plot(U_.'/Ud, (w_.'/(2*pi)), 'LineWidth', 1.2);
xlabel("U/U_D");
ylabel("\omega/2\pi [Hz]");



% %% pk method
% 
% % Aerodynamic matrix
% AICfun = @(k_) AIC_builder(k_, x_p, y_p, N_panels, c_ref, lambda, y0, b, M_inf);
% 
% % Init vel. vector
% dU = 5;
% U_ = dU:dU:1.6*Ud;
% 
% % Zero tol
% tol = 1e6;
% max_iter = 100;
% 
% % Init vars
% m_ = nan(N,N,length(U_));
% w_ = zeros(N,length(U_));
% g_ = zeros(N,length(U_));
% Umin = [];
% Umax = [];
% 
% %Initial guess
% w_(:,1) = fnat(i_nodes)*2*pi;
% 
% % Loop through vels
% for i = 1:length(U_)
% 
%     % Initilize timer
%     tic
% 
%     %Loop through modes
%     for j = 1:N
%         % Init conv vars
%         iter = 0;
%         res = tol;
% 
%         % Covergence loop
%         while iter < max_iter && res >= tol
%             % Update iter count
%             iter = iter +1;
% 
%             %Estim red freq
%             k = w_(j,i)*c/(2*U_(i));
% 
%             M_inf = U_(i)/a_inf;
% 
%             % Comput aerodin mat
%             Aefff = rho_inf/2*S_red*(AIC(k,M_inf)\(Ix_red + 1i*2*k/c*It_red));
% 
%             Keff = K_red - U_(i)^2*Aeff;
% 
%             % Eigs
%             I = eye(size(K_red));
%             O = zeros(size(K_red));
%             A = [Keff,O;O,I];
%             B = [O, -M_red, I, O];
%             [X,P] = eig(A,B);
%             p = diag(P);
% 
%             % Closes elem
%             [res,kast] = min(abs(real(p) - w_(j,i)*g_(j,i)) + abs(imag(p) - w_(j,i)));
% 
%             % Update the conv val
%             w_(j,i) = imag(p(kast(1)));
%             g_(j,i) = imag(p(kast(1)))/imag(p(kast(1)));
%             m_(:,j,i) = X(1:N,kast(1));
% 
% 
% 
%         end
% 
%         % Print
%         fprtintf("          Mode %i converged with %i iters",j,iter);
% 
% 
%         % Check stab
%         if i> 1 && g_(j,i-1) < tol && g_(j,i) > tol && isempty(Umin)
%             Umin = U_(i-1);
%             Umax = U_(i);
%         end
%     end
% 
%     % Update intial guess
%     if i<length(U_)
%         w_(:,i+1) = w_(:,i);
%         g_(:,i+1) = g_(:,1);
%     end
% 
% end
% 
% % Plots
% 
%     figure(3)
%     subplot(2,1,1)
%     hold on; box on;
%     plot(U_'/Ud,g_'*w_c);
%     xlabel("U/U_D");
%     ylabel("g");
% 
%     subplot(2,1,2)
%     hold on; box on;
%     plot(U_'/Ud,w_'/(2*pi));
%     xlabel("U/U_D");
%     ylabel("w/2pi");





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PROJECTE 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Now we can compute the pressure difference coefficients complex amplitudes
C_p = -(AIC\Wref);

% Pressure difference distribution.
% Com considerem un cas estàtic, el C_p actual ja és el C_p estàtic.

x_ = zeros(4, N_panels);   % 4 vèrtexs de cada panell
y_ = zeros(4, N_panels);
dy = b / Ny;

%Per visualitzar els panells amb Cp
for i = 1:Nx
    for j = 1:Ny
        k = Ny*(i-1) + j;

        % y-coordinates
        y1 = (j-1)*dy;
        y2 = j*dy;
        yc = (y1 + y2)/2;

        % local chord
        c1 = c_root * (1 - (1-lambda)*(y1/b));
        c2 = c_root * (1 - (1-lambda)*(y2/b));
        cc = c_root * (1 - (1-lambda)*(yc/b));

        % panel coordinates
        x_(1,k) = (i-1) * c1/Nx;
        y_(1,k) = y1;
        x_(2,k) =  i*c1/Nx;
        y_(2,k) = y1;
        x_(3,k) = i*c2/Nx;
        y_(3,k) = y2;
        x_(4,k) = (i-1)*c2/Nx;
        y_(4,k) = y2;
    end
end

% C_p físic per plotejar: part real
Cp_real = real(C_p);          % N_panels x 1

% Amplitud (per definir CLim de forma robusta)
Cp_amp = abs(C_p);
cpmax  = max(Cp_amp);

% Definim la superfície de l'ala (z=0)
%Zsurf = zeros(size(x_));
%figure("Units","normalized","Position",[0,0,0.5,1]);
%hold on; axis equal;

% Patch gris de l'ala (superfície neutra)
%p1 = patch(x_, y_, Zsurf,'FaceColor', 0.8*[1 1 1], 'EdgeColor', 'none');

% C_p per a cada panell com a fila 1xN_panels
Cp_plot = Cp_real.';   % transposat: 1 x N_panels

% Patch acolorit amb Cp
%p2 = patch(x_, y_, Zsurf,repmat(Cp_plot, 4, 1),'FaceColor', 'flat', 'EdgeColor', 'none');
%cb = colorbar;
%cb.Label.String = 'C_p';

%xlabel('x');
%ylabel('y');
%zlabel('z');
%title('Pressure difference distribution (static \alpha)');

%clim([-1, 1]*cpmax);
%view(12, 15);   % vista 3D; view(2) per planta
%grid on;



%%  STEP 6 - Calculations

% Lift i CL (CAS ESTÀTIC)
q_inf = 0.5 * rho_inf * U_inf^2; %p dinamica
S_panel = zeros(N_panels,1); % Àrea de cada panell 
for i = 1:Nx
    for j = 1:Ny
        k = Ny*(i-1) + j;
        S_panel(k) = (c_loc(j)/ Nx) * dy; %area del panell
    end
end
L_panel = Cp_real.* q_inf.* S_panel;  % Lift de cada panell
Lmat = reshape(L_panel, Ny, Nx);   % cada columna = mateixa x, cada fila = mateixa y
L_span = sum(Lmat, 2);            % sumem al llarg de la corda (direcció x)
L  = 2*sum(L_panel); %Lift DE TOT EL PLANEJADOR
CL = L / (q_inf * S); 

fprintf('\n--- Resultats aerodinàmics ---\n');
fprintf('Lift total  L  = %.3f N\n', L);
fprintf('Coeficient CL  = %.4f\n', CL);

%------------------------------------------
% Massa total i Mass Ratio(mu)

y_span = span.y; % span.y va de 0 a b/2
m_span = span.m; % span.m és massa per unitat de longitud [kg/m]

M_semi = trapz(y_span, m_span); % Massa semiala (integració numerica)
M_tot = 2 * M_semi;    % Massa ala completa [kg]

% Mass ratio mu (ENUNCIAT)
mu = M_tot/(pi*rho_inf*c_root*S);

fprintf('\n--- Massa ---\n');
fprintf('Massa total ala     = %.3f kg\n', M_tot);
fprintf('Mass ratio mu       = %.4f\n', mu);

%------------------------------------------
% Lift-to-mass ratio

fprintf('\n--- Lift-to-mass ratio ---\n');
fprintf('CL / mu = %.4f\n', CL/mu);

%------------------------------------------
% Stress and Stress-to-mass Ratio

E = material.YoungModulus;
G  = E/(2*(1+material.Poisson));
y = span.y;
c_y = span.c_y;
Ny = length(y);

% Recuperar el primer mode (1ra frequencia)
Phi = modes.Phi; % modes propis (només DOFs lliures)
ndof = max(modes.freeDOFs); %n total de DOFs
mode1 = Phi(:,1); % primer mode            
mode_full = zeros(ndof,1);
mode_full(modes.freeDOFs) = mode1;

zeta_y  = zeros(Ny,1);
theta_y = zeros(Ny,1);

for i = 1:Ny
    k = (i-1)*3;
    zeta_y(i)  = mode_full(k + 2); %DOF 2 ->zeta
    theta_y(i) = mode_full(k + 3); %DOF 3 ->theta
end

% Derivades en y
dzeta_dy  = gradient(zeta_y, y);
dtheta_dy = gradient(theta_y, y);

% h_max i r_sc_max escalats amb la corda
t_root = param(2); % thickness root
t_y = t_root * (c_y / c_y(1)); % thickness(y)
h_max = t_y / 2; % fibra extrema

xsc_root = section.xsc(1);  % shear center root
r_sc_max = abs(xsc_root) * (c_y / c_y(1)); % distancia maxima SC(y)

% Tensió (SEGONS ENUNCIAT)
sigma_y = sqrt((E/2.* h_max.* dzeta_dy).^2 +(3*G.* r_sc_max.* dtheta_dy).^2 );
sigma_max = max(sigma_y);

% Coeficient C_sigma i ratio final
q_inf = 0.5*rho_inf*U_inf^2;
C_sigma = sigma_max/(q_inf*span.S);
Csigma_over_mu = C_sigma / mu;

fprintf('\n--- Stress ---\n');
fprintf('Sigma max        = %.3f MPa\n', sigma_max/1e6);
fprintf('C_sigma          = %.4f\n', C_sigma);
fprintf('C_sigma/mu       = %.4f\n', Csigma_over_mu);


% -----------------------
%Guardem els resultats del cas

metrics.y_panel = y_mid;
metrics.b = b;
metrics.lambda = lambda;
metrics.Nx = Nx;
metrics.Ny = Ny;
metrics.S = S;
metrics.AR = AR;
metrics.Cp = Cp_plot;
metrics.f1 = f1;
metrics.L = L;
metrics.L_span = L_span;
metrics.CL = CL;
metrics.mass = M_tot;
metrics.mu = mu;
metrics.CL_over_mu = CL/mu;
metrics.sigma_max = sigma_max;
metrics.C_sigma = C_sigma;
metrics.Csigma_over_mu = Csigma_over_mu;

metrics.span     = span;
metrics.modes    = modes;
metrics.c_root   = c_root;
metrics.x_panel  = x_panel;
metrics.y_panel  = y_panel;
metrics.material = material;


end
