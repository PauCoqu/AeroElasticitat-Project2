%% CODES COPIED FROM OTHER CLASSES FALTEN COSES AQUI DALT, SUPOSO QUE DEL CODI DE LA CLASSE 3

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


% BOundary conditions
I_fix = [
    1:3*(Ny+1):N_dof,
    2:3*(Ny+1):N_dof,
    3:3*(Ny+1):N_dof,
    ];

I_free = setdiff(1:N_dof,I_fix);

%Structural modes
q_mod = zeros(N_dof,length(I_free));
[q_mod(I_free,:),w2] = eig(K(I_free,I_free),M(I_free,I_free));
fnat = sqrt(diag(w2)/(2*pi));

% Model reduction
i_modes = [1,2,3,5];

N = length(i_modes);

% Reduced matrices
M_red = q_mod(:,i_modes)'*M*q_mod(:,i_modes);
K_red = q_mod(:,i_modes)'*K*q_mod(:,i_modes);
S_red = q_mod(:,i_modes)'*S;
It_red = It*q_mod(:,i_modes);
Ix_red = Ix*q_mod(:,i_modes);


%Flow properties
rho_inf = 1.25;
a_inf = 343;
Ud = 145; %Divergence speed

%% p method

% Aerodynamic matrices
k = 0; % % quasi-steady

M_inf = 0; % Incompressibility

AIC = dlmAIC_vec(x_p,y_p,k,M_inf,c/2);

A0_red = rho_inf/2*S_red*(AIC\Ix_red); % Parenthesis important!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
A1_red = -rho_inf/2*S_red*(AIC\It_red);


% Initialize velocity vector
dU = 5;
U_ = dU:dU:1.6*Ud;

% Zero tolerance
tol = 1e-6;

% Initialize variables
p_ = nan(2*N,length(U_));
m_ = nan(N,2*N,length(U_));
U_min = [];
U_max = [];

% Loop through velocities
for i = length(U_)

    % Effecftive matrices
    Keff = K_ref - U_(i)^2*A0_red;
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
        m_(:,:,i) = X(1:N,:);
    else
        m2sort = 1:2*N;
        for j = 1:length(p)
            [~,jmin] = min(abs(real(p_(j,i-1))-real(p(m2sort)))+ abs(imag(p_(j,i-1))-imag(p(m2sort))));

            p_(j,i) = p(m2sort(jmin));
            m_(:,j,i) = X(1:N,m2sort(jmin));
            m2sort = setdiff(m2sort,m2sort(jmin));

        end
    end

    %Check stability
    if max(real(p_(:,i))) > tol && isempty(Umin) % Check if the min vel is encounterd. If yes, the flutter speed is betweem this iteration and the past one
        Umin = U_(i-1);
        Umax = U_(i);
        break;
    end


    %% PLOTS

    figure
    subplot(2,1,1)
    hold on; box on;
    plot(U_/Ud,real(p_)*c./(2*U));
    xlabel("U/U_D");
    ylabel("p_Rc/2U");

    subplot(2,1,2)
    hold on; box on;
    plot(U_/Ud,imag(p_)*c./(2*pi));
    xlabel("U/U_D");
    ylabel("p_I/2\pi");






%% k method
% We need to recompute the matrices at each iteration

% Aerodynamic
M_inf = 0;
AIC = @(k_) dlmAIC_vec(x_p,y_p,k_,M_inf,c/2);

% Initialized reduced freq
dinvk = 0.1;
invk_ = dinvk:dinvk:20;

% Zero tolerance
tol = 1e-6;

% Initialize variables
m_ = nan(N,N,length(invk_));
l_ = nan(N,length(invk_));
Umin = [];
Umax = [];

% Loop through reduced frequencues
for i = 1:length(invk_)

    % Initilize timer
    tic



    % Effective matrices
    k = 1/invk-(i);
    Beff = M_red + rho_inf/2*S_red*(AIC(k)\(1i*c/(2*k)*It_red + c^2/(4*k^2)*Ix_red));

    %Eigenvalues
    [X, L] = eig(K_red,Beff);
    l = 1./diag(L);


    % Sort modes
    if i == 1
        l_(:,i) = l;
        m_(:,:,i) = X;
    else
        m2sort = 1:N;
        for j = 1:length(l)
            [~,jast] = min(abs(real(l_(j,i-1))-real(l(m2sort))) + abs(imag(l_(j,i-1))-imag(l(m2sort))));
            l_(j,i) = l(m2sort(jast));
            m_(:,j,i) = X(:,m2sort(jast));
            m2sort = setdiff(m2sort,m2sort(jast));
        end
    end

    %Check stability
    if max(imag(l_(:,i))) > tol && isempty(Umin)

        [~,jast] = max(imag(l_(:,i)));
        Umin = sqrt(1./real(l_(jmax,i-1)))*c(2*invk_(i-1);
        Umax = sqrt(1./real(l_(jmax,i)))*c(2*invk_(i-1);
    end

    % Print iteration time
    %ja tenim el toc, ha jo fare aixo
end


% Recover parameters
g_ = imag(l_)./real(l_);
w_ = sqrt(1./real(l_));
U_ = w_*c/2.*invk_;

% Plots

    figure(2)
    subplot(2,1,1)
    hold on; box on;
    plot(U_'/Ud,g_');
    xlabel("U/U_D");
    ylabel("g");

    subplot(2,1,2)
    hold on; box on;
    plot(U_'/Ud,w_'/(2*pi));
    xlabel("U/U_D");
    ylabel("w/2pi");


%% pk method

% Aerodynamic matrix
AIC = @(k_,M_) dlmAIC_vec(x_p,y_p,k_,M_inf,c/2);

% Init vel. vector
dU = 5;
U_ = dU:dU:1.6*Ud;

% Zero tol
tol = 1e6;
max_iter = 100;

% Init vars
m_ = nan(N,N,length(U_));
w_ = zeros(N,length(U_));
g_ = zeros(N,length(U_));
Umin = [];
Umax = [];

%Initial guess
w_(:,1) = fnat(i_nodes)*2*pi;

% Loop through vels
for i = 1:length(U_)

    % Initilize timer
    tic

    %Loop through modes
    for j = 1:N
        % Init conv vars
        iter = 0;
        res = tol;

        % Covergence loop
        while iter < max_iter && res >= tol
            % Update iter count
            iter = iter +1;
            
            %Estim red freq
            k = w_(j,i)*c/(2*U_(i));

            M_inf = U_(i)/a_inf;

            % Comput aerodin mat
            Aefff = rho_inf/2*S_red*(AIC(k,M_inf)\(Ix_red + 1i*2*k/c*It_red));

            Keff = K_red - U_(i)^2*Aeff;

            % Eigs
            I = eye(size(K_red));
            O = zeros(size(K_red));
            A = [Keff,O;O,I];
            B = [O, -M_red, I, O];
            [X,P] = eig(A,B);
            p = diag(P);

            % Closes elem
            [res,kast] = min(abs(real(p) - w_(j,i)*g_(j,i)) + abs(imag(p) - w_(j,i)));

            % Update the conv val
            w_(j,i) = imag(p(kast(1)));
            g_(j,i) = imag(p(kast(1)))/imag(p(kast(1)));
            m_(:,j,i) = X(1:N,kast(1));

        

        end

        % Print
        fprtintf("          Mode %i converged with %i iters",j,iter);


        % Check stab
        if i> 1 && g_(j,i-1) < tol && g_(j,i) > tol && isempty(Umin)
            Umin = U_(i-1);
            Umax = U_(i);
        end
    end

    % Update intial guess
    if i<length(U_)
        w_(:,i+1) = w_(:,i);
        g_(:,i+1) = g_(:,1);
    end

end

% Plots

    figure(3)
    subplot(2,1,1)
    hold on; box on;
    plot(U_'/Ud,g_'*w_c);
    xlabel("U/U_D");
    ylabel("g");

    subplot(2,1,2)
    hold on; box on;
    plot(U_'/Ud,w_'/(2*pi));
    xlabel("U/U_D");
    ylabel("w/2pi");

