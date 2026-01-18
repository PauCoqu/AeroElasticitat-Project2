function metrics = wing_case(geometry, param, density,youngModulus,poissonRatio, b, lambda, Nx, Ny,U_inf,y0, i_modes,alpha_deg)
%% STEP 1: Section properties (del pefil) ('beam_properties')

% Section properties
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
    
ctip = lambda * c_root;    % corda punta
Surf = b * (c_root + ctip); % Area de les DOS ales
AR   = (2*b)^2 / Surf;      % Aspect ratio 

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
alpha = alpha_deg*pi/180;  %rad

% Reduced frequency -> Frequencia d'oscilació (k>0 -> flutter i altres, k=0 -> CAS ESTÀTIC)
k = 0; %2*pi*freq*c/(2*U_inf)
M_inf = 0; % Mach = U_inf/a_inf; %0; %CAS INCOMPRESSIBLE;

% Wref = alpha*ones(N_panels,1); PROJECTE 1
%Wref = -1i*k*eta/(c/2)*ones(N_panels,1); 

% AIC matrix coefficients (com buildAIC.m)
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
% % Ud -> 1a inestabilitat -> det(K - U^2*A0)=0 -> K*phi = U^2*A0*phi -> Ud = sqrt(menor U^2 real i positiu)


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

A0 = 0.5*rho_inf * (S*(AIC\Ix));
A0_free = A0(I_free,I_free);      % fins aquí hauria d'estar tot bé

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute divergence speed with ALL MODES
[X_div,U_2inf] = eigs(sparse(K(I_free,I_free)), sparse(A0(I_free,I_free)), 10, 'smallestabs');

autov = diag(U_2inf);  %extra         % lambda(autov) = U^2
tol = 1e-8 * max(1, max(abs(autov)));
valid = abs(imag(autov)) < tol & real(autov) > 0;
U_inf2 = sort(sqrt(real(autov(valid))));
U_inf_verif = sort(diag(sqrt(U_2inf)))
if isempty(U_inf2)
    disp('No divergence speed found');
end
Ud_allmodes = U_inf2(1) %DIVERGENCE SPEED, valor real


% Compute divergence speed with REDUCED SYSTEM (i_modes)

q_mod = zeros(N_dof,length(I_free));
[q_mod(I_free,:),w2] = eig(K_free,M_free);
freq = sqrt(diag(w2))/(2*pi); %frequencia 
K_red = q_mod(:,i_modes)'*K*q_mod(:,i_modes);
A0_red = q_mod(:,i_modes)'*A0*q_mod(:,i_modes);

[X_div,U_3inf] = eig(K_red,A0_red);

lambda2 = diag(U_3inf);  %extra         % lambda = U^2
tol = 1e-8 * max(1, max(abs(lambda2)));
valid = abs(imag(lambda2)) < tol & real(lambda2) > 0;
U_inf3 = sort(sqrt(real(lambda2(valid))));

%U_inf3 = sort(diag(sqrt(U_3inf)));

%Ud_reduced = U_inf3(1) %DIVERGENCE SPEED, assumint que el 1r eigenvalue es real

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ud = Ud_allmodes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 2. RISK OF FLUTTER (ara mateix ens depen de Ud. Sino tenim Ud peta)
k = 0; % quasi-steady  
M_inf = 0; % Incompressibility %JA DEFINIT MÉS AMUNT

%Reuse modal results from modal_analysis_wing
q_mod = zeros(N_dof, size(modes.Phi_free,2));
q_mod(I_free,:) = Phi_free;

% Select modes for p-method as we do a Model reduction
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
U_ = 0.5:0.5:1.6*Ud;

% Initialize variables
p_ = nan(2*N_modes,length(U_));
k_post = nan(size(U_));
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
    if i > 1
    r_prev = max(real(p_(:,i-1)));
    if isempty(U_min) && (r_prev <= 0) && (max(real(p_(:,i))) > 0)
        disp(' Creuem el 0 entre U_min i U_max ');
        U_min = U_(i-1)
        U_max = U_(i)
    end
    end
    
  
%%%%%%%%%%%%%%%%%%%%%%
%Check validity p method (k<<)
tol_w = 1e-6;                         
osc = abs(imag(p_(:,i))) > tol_w;
if any(osc)
    [~, jj] = max(real(p_(osc,i)));  % most unstable oscillatory root
    idx = find(osc);
   jcrit = idx(jj);
    omega = abs(imag(p_(jcrit,i)));    % rad/s
    k_validation(i) = omega * c_root / (2*U_(i)); %reduced frequency
end
%%%%%%%%%%%%%%%%%%%%%%


end

%% PLOTS

% --- Subplot real parts ---
figure 
subplot(2,1,1) 
hold on; box on;
plot(U_/Ud, real(p_).*c_root./(2*U_), 'LineWidth', 2); yline(0,'k--')
xlabel("U_{\infty}/U_D"); ylabel("p_R c / 2U");
subplot(2,1,2)
hold on; box on;
plot(U_/Ud, imag(p_)/(2*pi), 'LineWidth', 2); yline(0,'k--')
xlabel("U_{\infty}/U_D"); 
ylabel("p_I / 2\pi");

% Validació de k (p method)
figure; hold on; box on;
plot(U_, k_validation, 'LineWidth', 2);
yline(0.05,'k--'); 
yline(0.1,'k--');
xlabel('U_{\infty} [m/s]');
ylabel('k = \omega c/(2U)');
%xlim([90 120])   %<-- limit U range
%ylim([0 10])   %<-- limit U range




%--------------------------
%% GUSTS    

% Gust data (enunciat)
sigma_g = 1; % [m/s]
k0 = 10^(-28/9);

% Speed sweep for gust analysis
U_g = 5:15:Ud*1.4; % [m/s] up to divergence

% Frequency vector (igual q a classe, logspace per low freqs, then linear)
f = [0, logspace(-5,0,50), 1.2:0.2:20];
w = 2*pi*f;

% C_red = 0
C_red = zeros(size(M_red));


% Tip vertical acceleration at tip shear center

% Tip index
x_sc_line = span.x_sc(1); % shear center x at root
y_tip = y(end);
tip_nodes = find(abs(malla.nodes(:,2) - y_tip) < 1e-12);
[~,ii] = min(abs(malla.nodes(tip_nodes,1) - x_sc_line));
i_tip = tip_nodes(ii);

% Acceleration measurement vector
Ie = zeros(N_dof,1);
Ie(3*(i_tip-1) + 1) = 1; % eta = 1 al tip node
Ie_red = Phi.' * Ie;

% Allocate outputs
Phi_gg = zeros(length(f), length(U_g));
Phi_aa = zeros(length(f), length(U_g));
Phi_ss = zeros(length(f), length(U_g));
H_acc = zeros(length(f), length(U_g));
H_sig = zeros(length(f), length(U_g));

for iU = 1:length(U_g)
    U_inf_i = U_g(iU);
    for jf = 1:length(f)
        if f(jf) == 0, continue; end  % skip f=0 per evitar divs per 0
        omega = w(jf);
        k = omega * c_root / (2 * U_inf_i);

        % AIC at this k unsteady (la del profe)
        AICk = buildAIC(x_p, y_p, c_root, lambda, y0, b, M_inf, k);

        %Aeroelastic coupling matrices
        Aeff_full = rho_inf/2 * ( S * (AICk \ ( 1i*k*2/c_root * It + Ix )) );

        % Gust loads vector
        Ig = exp(1i*k*2/c_root * x_p(:,3));
        Ag_full = -rho_inf/2 * ( S * (AICk \ Ig) );

        % Reduce to modal coordinates
        A_red = Phi.'* Aeff_full * Phi;
        Ag_red = Phi.'* Ag_full;
        
        % System FRF: modal response per unit gust velocity amplitude Wg
        % Q = (-w^2 M + i w C + K - U^2 A_red)^(-1) * Ag_red
        Den = (-omega^2)*M_red + 1i*omega*C_red + K_red - (U_inf_i^2)*A_red;
        Q_perWg = Den \ Ag_red;

        % tip vertical acceleration FRF
        eta_tip_perWg = Ie_red.'* Q_perWg;
        H_acc(jf,iU) = (-omega^2) * eta_tip_perWg;

        % sigma_max FRF Reconstruct full DOFs: q_full = Phi * Q
        q_full_perWg = Phi * Q_perWg;
        sigma_max_perWg = sigmaMax_fromResponse(q_full_perWg, malla.nodes, span, material, x_sc_line);
        H_sig(jf,iU) = sigma_max_perWg;

        % Gust PSD in k
        Phi_k = gustPSD_k(k, sigma_g, k0);
        Phi_gg(jf,iU) = Phi_k * (c_root / (2 * U_inf_i));

        % Response PSDs
        Phi_aa(jf,iU) = abs(H_acc(jf,iU))^2 * Phi_gg(jf,iU);
        Phi_ss(jf,iU) = abs(H_sig(jf,iU))^2 * Phi_gg(jf,iU); % diapo 199
    end
    fprintf(' Gust: %03i/%03i (U=%.1f m/s)\n', iU, length(U_g), U_inf_i);
end

% RMS values
g_RMS = sqrt(trapz(w, Phi_gg,1)); % RMS gust velocity
a_RMS = sqrt(trapz(w, Phi_aa,1)); % RMS tip acceleration
s_RMS = sqrt(trapz(w, Phi_ss,1)); % RMS sigma_max

% Acceleration PSD curve (picked a few speeds)
U_plot = U_g(round(linspace(1, length(U_g), 4)));
[~,idxU] = ismember(U_plot, U_g);
figure; hold on; box on;
for kU = idxU
    plot(f, Phi_aa(:,kU), 'LineWidth', 1.5);
end
set(gca,'XScale','log','YScale','log');
xlabel('f [Hz]'); ylabel('\Phi_{\eta} [(m/s^2)^2/(rad/s)]');
title('Tip acceleration PSD (selected U_\infty)');
legend(arrayfun(@(U) sprintf('U=%.0f m/s',U), U_plot, 'UniformOutput', false), 'Location','best');

% Stress PSD curve
figure; hold on; box on;
for kU = idxU
    plot(f, Phi_ss(:,kU), 'LineWidth', 1.5);
end
set(gca,'XScale','log','YScale','log');
xlabel('f [Hz]'); ylabel('\Phi_{\sigma\sigma} [Pa^2/(rad/s)]');
title('Max stress PSD (selected U_\infty)');
legend(arrayfun(@(U) sprintf('U=%.0f m/s',U), U_plot, 'UniformOutput', false), 'Location','best');

% Acceleration RMS vs speed
figure; hold on; box on;
plot(U_g, a_RMS, 'LineWidth', 2);
xlabel('U_\infty [m/s]'); ylabel('\eta_{RMS} [m/s^2]'); % Hauria de ser \ddot{\eta}_{RMS} [m/s^2] pero no sé posar ddot
title('Tip acceleration RMS vs U_\infty');

%Stress RMS vs speed
figure; hold on; box on;
plot(U_g, s_RMS, 'LineWidth', 2);
xlabel('U_\infty [m/s]'); ylabel('\sigma_{max,RMS} [Pa]');
title('Max stress RMS vs U_\infty');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PROJECTE 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analisi aerodinamic

%  Cp DLM (estatic)
alpha = deg2rad(5);            % rad
Wref  = alpha * ones(N_panels,1);
C_p   = -(AIC\Wref);
Cp_real = real(C_p);

x_ = zeros(4,N_panels);
y_ = zeros(4,N_panels);

for i = 1:Nx
    for j = 1:Ny
        k  = Ny*(i-1) + j;
        y1 = y(j);
        y2 = y(j+1);
        c1 = c_root * (1 - (1-lambda) * (y1 - y0)/b);
        c2 = c_root * (1 - (1-lambda) * (y2 - y0)/b);
        x_(1,k) = (i-1)*c1/Nx;  y_(1,k) = y1;
        x_(2,k) =  i   *c1/Nx;  y_(2,k) = y1;
        x_(3,k) =  i   *c2/Nx;  y_(3,k) = y2;
        x_(4,k) = (i-1)*c2/Nx;  y_(4,k) = y2;
    end
end

% --- Lift and CL ---
q_inf = 0.5*rho_inf*U_inf^2;

S_panel = zeros(N_panels,1);
for k = 1:N_panels
    dy_loc = y_(4,k) - y_(1,k);
    dx1 = x_(2,k) - x_(1,k);
    dx2 = x_(3,k) - x_(4,k);
    S_panel(k) = 0.5*(dx1+dx2) * dy_loc;     % area trapezoidal
end

L_panel = Cp_real .* q_inf .* S_panel;
L_semi  = sum(L_panel);                     % semiwings lift
L_total = 2*L_semi;                         % total lift (two wings)
ctip    = lambda*c_root;
S_total = b*(c_root + ctip);                % total area (2 wings)
CL      = L_total / (q_inf * S_total);

fprintf('\n--- Resultats aerodinàmics ---\n');
fprintf('Lift semiala   = %.3f N\n', L_semi);
fprintf('Lift total     = %.3f N\n', L_total);
fprintf('CL             = %.4f\n', CL);
fprintf('max|Cp|        = %.3f\n', max(abs(Cp_real)));


%% COMPARACIÓ x_cp amb x_sc

x_cent = mean(x_(1:4,:), 1).';
L_sec    = zeros(Ny,1);
M_LE_sec = zeros(Ny,1);
y_mid    = zeros(Ny,1);
c_mid    = zeros(Ny,1);

for j = 1:Ny
    idx = j:Ny:N_panels;

    L_sec(j)    = sum(L_panel(idx));
    M_LE_sec(j) = sum(L_panel(idx) .* x_cent(idx));

    y_mid(j) = 0.5*( y_(1,idx(1)) + y_(4,idx(1)) );
    c_mid(j) = c_root * (1 - (1-lambda) * (y_mid(j) - y0)/b);
end

% Aerodynamic center per seccio (quarter-chord)
x_ac_sec = 0.25 * c_mid;

% x_cp (x_cp = M_LE/L)
x_cp_sec = M_LE_sec ./ L_sec;
%comprovació (x_cp = xac - M_ac/L)
M_ac_sec   = M_LE_sec - x_ac_sec .* L_sec;
x_cp_check = x_ac_sec + M_ac_sec ./ L_sec;
fprintf('max|x_cp - x_cp_check| = %.3e\n', max(abs(x_cp_sec - x_cp_check)));

% x_sc
x_sc_sec = span.x_sc(1:Ny);

figure; hold on; box on;
plot(y_mid, x_cp_sec, 'LineWidth', 2);
plot(y_mid, x_sc_sec, 'LineWidth', 2);
xlabel('y [m]'); ylabel('x [m]');
legend('x_{cp}','x_{sc}','Location','best');
title('Comparacio x_{cp} vs x_{sc} (metres, LE=0)');

%% ALTRES PARAMETRES AERODINAMICS

% -------------------------------
% Massa total i Mass Ratio(mu)

y_span = span.y; % span.y va de 0 a b/2
m_span = span.m; % span.m és massa per unitat de longitud [kg/m]

M_semi = trapz(y_span, m_span); % Massa semiala (integració numerica)
M_tot = 2 * M_semi;    % Massa ala completa [kg]

% Mass ratio mu (ENUNCIAT)
mu = M_tot/(pi*rho_inf*c_root*Surf);

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

metrics.CL_over_mu = CL/mu;
metrics.sigma_max = sigma_max;
metrics.C_sigma = C_sigma;
metrics.Csigma_over_mu = Csigma_over_mu;

end