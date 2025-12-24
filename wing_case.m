function metrics = wing_case(density,youngModulus,poissonRatio, b, lambda, Nspan, Nx, Ny, alpha_deg,U_inf)

%% STEP 1: Section properties (del pefil) ('beam_properties')  [section]

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
S    = b * (c_root + ctip); % Area de les DOS ales
AR   = (2*b)^2 / S; % Aspect ratio 

% Corda a cada punt: c(y) = c_root*[1 - (1-lambda)*2y/b]
y     = linspace(0, b, Ny)'; % coordenades de cada punt
c_y = c_root * (1 - (1 - lambda) * y/b); %corda cada perfil

%Calcul de les propietats al llarg de l'envergadura [span]
[span] = spanwise_distributions(section, material, c_y, c_root, y, b, S, AR, lambda);

%% STEP 3 - Euler-Bernoulli beam 

% AQUEST PAS NOMÉS CAL EXECUTAR-LO UNA VEGADA!!!
% En la següent funció (beamElementSymbolic):
%S'obtenen les matrius: M_e (masa) ,K_e (rigidesa),B_e (cargues aerodinamiques)
% MATLAB genera les funcions beamMass.m, beamStiffness.m y beamForce.m.
%[M_e, K_e, B_e] = beamElementSymbolic();

% ARA JA TENIM ELS ARXIUS BEAM STIFFNESS, BEAMMASS, BEAMFORCE

%% STEP 4 - Modal Analysis 

[modes] = modal_analysis_wing(span, material); %Calculem les matrius
f1 = modes.f(1);   % primera freq pròpia

%% STEP 5 - DOUBLET LATTICE 
% Ja podem correr el solver | %run('doublet_lattice.mlx');

% Discretització superfície aerodinàmica (model DLM 2D x-y).
N_panels = Nx*Ny; %totals

% Geometry data que falta (pels panells aerodinamics)
y_mid = ((1:Ny)-0.5)*b/Ny;   % Centre del panell
c_loc = c_root * (1-(1-lambda)*y_mid/b);

% Let's define the coordinates of each relevant point in all panels
x_panel = zeros(N_panels,3);
y_panel = zeros(N_panels,3);
dy = b/Ny;

for i = 1:Nx
    for j = 1:Ny
        k = Ny*(i-1)+j;

        % Y-coordinates
        y1 = (j-1)*dy; %costat inferior
        y2 = j*dy; %costat superior
        yc = (y1 + y2)/2; %punt mig

        % Chord
        c1 = c_root * (1 - (1-lambda)*(y1/b)); %costat inferior
        c2 = c_root * (1 - (1-lambda)*(y2/b)); %costat superior
        cc = c_root * (1 - (1-lambda)*(yc/b)); %punt mig

        % DOUBLET segment endpoints
        x1 = (i-3/4) * c1 / Nx; 
        x2 = (i-3/4) * c2 / Nx;
        xc = (i-1/4) * cc / Nx;

        x_panel(k,1) = x1;    
        y_panel(k,1) = y1;
        x_panel(k,2) = x2;  
        y_panel(k,2) = y2;
        x_panel(k,3) = xc;   
        y_panel(k,3) = yc;
    end
end


% Visualize mesh
%figure("Units","normalized","Position",[0,0,1,1]);
%hold on; axis equal;
%patch([0, 0, ctip, c_root],[0,b,b,0], zeros(1,4), 'FaceColor','none', 'EdgeColor','k');
% Línies dels panells (vortex rings)
%patch(x_panel(:,1:2)', y_panel(:,1:2)', zeros(2,N_panels),'FaceColor','none','EdgeColor','b');
%plot(x_panel(:,3), y_panel(:,3), '.r');
%xlabel x; ylabel y;


%% STEP 5.1 - Aerodynamic influence coefficients

% Fluid data
rho_inf = 1.25;
%a_inf = 343;
%eta = 0.1*(c/2);  %Amplitud del moviment vertical (CAS ESTATIC)
%freq = 1; % Hz (CAS ESTATIC)
alpha = alpha_deg*pi/180;  %rad

% Reduced frequency -> Frequencia d'oscilació (k>0 -> flutter i altres, k=0 -> CAS ESTÀTIC)
k = 0; %2*pi*freq*c/(2*U_inf);
% Mach
M = 0; %U_inf/a_inf; %0; %CAS INCOMPRESSIBLE;


% Wref = alpha*ones(N_panels,1); PROJECTE 1
%Wref = -1i*k*eta/(c/2)*ones(N_panels,1); 

%PROJECTE 2--------------------------



%------------------------------------
% AIC matrix coefficients
AIC = zeros(N_panels,N_panels);
for i = 1:N_panels
    % Collocation point coordinates of panel "i"
    x_i = x_panel(i,3);
    y_i = y_panel(i,3);
    % Loop through all doublet segments
    for j = 1:N_panels
        % We add the induced velocity contribution to the AIC
        AIC(i,j) = AIC(i,j) + w_doublet(x_i,y_i,x_panel(j,:),y_panel(j,:),c_root/2,M,k,1);
        % If an image panel is present
        % 1) The symmetry makes the orientation reversed (point 1->2 and 2->1)
        % 2) The symmetry makes span y-coordinate negative
        AIC(i,j) = AIC(i,j) + w_doublet(x_i,y_i,x_panel(j,[2,1,3]),-y_panel(j,[2,1,3]),c_root/2,M,k,1);
    end
end

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
metrics.Nspan = Nspan;
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
