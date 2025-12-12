clear all;
close all;

c = 0.75;
b = 2;
h = 0.01;

E = 70e9;
nu = 0.3;
rho = 2300;

%Number of elements
Nx = 20;
Ny = 40;

% Structural mesh
nodes = zeros((Nx+1)*(Ny+1),2);

% Coordinates of the nodes
for i = 1:Nx+1
    for j = 1:Ny+1
        k = (Ny+1)*(i-1)+j;
        nodes(k,1) = (i-1)*c/Nx;
        nodes(k,2) = (j-1)*h/Ny;
    end
end

% Element nodes
elem = zeros(Nx*Ny,4);
for i = 1:Nx
    for j = 1:Ny
        e = Ny*(i-1) + j;
        elem(e,1) = (Ny+1)*(i-1) + j;
        elem(e,2) = (Ny+1)*i + j;
        elem(e,3) = (Ny+1)*i + j+1;
        elem(e,4) = (Ny+1)*(i-1) + j+1;
    end
end


%Aerodynamic mesh
x_p = zeros((Nx+1)*Ny,5);
y_p = zeros((Nx+1)*Ny,5);

for i = 1:Nx+1
    for j = 1:Ny+1
        k = Ny*(i-1) + j;

        %Corner 1
        x_p(k,1) = (i-3/4)*c/Nx;
        y_p(k,1) = (j-1)*b/Ny;
        %Corner 2
        x_p(k,2) = (i-3/4)*c/Nx;
        y_p(k,3) = j*b/Ny;
        %Corner 3
        if i>Nx
            x_p(k,3) = x_p(k,2) + 20*c;
        else
            x_p(k,3) = (i+1/4)*c/Nx;
        end
        y_p(k,1) = j*b/Ny;
        %Corner 4
        if i>Nx
            x_p(k,4) = x_p(k,1) + 20*c;
        else
            x_p(k,4) = (i+1/4)*c/Nx;
        end
        y_p(k,4) = (j-1)*b/Ny;

        if k>Nx
            x_p(k,5) = nan;
            y_p(k,5) = nan;
        else
            x_p(k,5) = (i-1/4)*c/Nx;
            y_p(k,5) = (j-1/2)*b/Ny;

        end
 
    end
end

%Aerodynamic influence coefficients
%[AIC,AIC_w] = vlmAIC(x_p,y_p,Nx,Ny,1);

%Kutta condition (steady case)
I_w = zeros(Ny,Ny*Nx);
for j = 1:Ny
    I_w(j,Ny*(Nx-1)+j) = 1;
end
%AIC = AIC + AIC_w*I_w;

%Assembly of structural matrices
N_dof = 3*size(nodes,1);
M = zeros(N_dof,N_dof);
K = zeros(N_dof,N_dof);
S = zeros(N_dof,Nx*Ny);

for i = 1:Nx*Ny
    % Element size
    a_e = (nodes(elem(i,2),1) - nodes(elem(i,1),1))/2;
    b_e = (nodes(elem(i,4),2) - nodes(elem(i,1),2))/2;
    % Element matrices
    M_e = plateMass(a_e,b_e,h,rho);
    K_e = plateStiffness(a_e,b_e,h,E,nu);
    S_e = plateForce(a_e,b_e);

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


%Structural
q_mod = zeros(N_dof,size(M_free,1));

[q_mod(I_free,:),w2] = eig(K_free,M_free);
freq = sqrt(diag(w2))/(2*pi);

%Plot first 6 nodes
%Foto al mobil IMG_5887.HEIC
x_ = nodes(:,1);
y_ = nodes(:,2);

figure;

for i = 1:6
    z_ = q_mod(1:3:N_dof,1);
    subplot(2,3,i);
    hold on;
    box on;
    patch(x_(elem'),y_(elem'),z_(elem'),zeros(size(elem')),...
        'FaceColor','flat','EdgeColor','k');
    xlabel('x');
    ylabel('y');
    zlabel('\eta');
    axis equal;
    axis tight;
    view(50,30);
    title(sprintf('f = %.3g Hz',freq(i)));
end


%Flow
rho_inf = 1.2;
U_inf = 20;
alpha = 15*pi/180;

%Velocity vector
W_ref = U_inf*alpha*ones(Nx*Ny,1);

%Pressure difference
D = zeros(Nx*Ny,Nx*Ny);
for i = 1:Nx
    for j = 1:Ny
        D(Ny*(i-1)+j,Ny*(i-1)+j) = 1/(c/Nx);
        if i>1
            D(Ny*(i-1)+j,Ny*(i-2)+j) = -1/(c/Nx);
        end

    end
end
delta_p = rho_inf*U_inf*D*(-AIC\W_ref);

%Deflections of the panel
q_aoa = zeros(N_dof,1);
q_aoa(I_free) = K_free\(S_free*delta_p);

%Plot results
figure;
z_ = q_aoa(1:3:N_dof);
   patch(x_(elem'),y_(elem'),z_(elem'),[1;1;1;1]*delta_p',...
        'FaceColor','flat','EdgeColor','k');
    xlabel('x');
    ylabel('y');
    zlabel('\eta');
    axis equal;
    axis tight;
    col = colorbar('south');
    col.Label.String = '\Deltap [N/m^(2)]';
    view(50,30);
    title("Deflection [m]");
