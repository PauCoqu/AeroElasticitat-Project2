%% Project 1 - Wing design and modelling
% Pau Cornudella Quer
% 220351 - ADVANCED AEROELASTICITY
%-------------------------------------------
clear; clc; close all

%% STEP 1: Section properties (del pefil) ('beam_properties')  [section]

%Material
density = 2300; %rho Alumini (ENUNCIAT)
youngModulus = 69000000000; % E alumini (ENUNCIAT)
poissonRatio = 0.3; % possion alumini (ENUNCIAT)


%% STEP 2 - Spanwise distributions

% VALORS DE DISSENY
b_vec = [12] ; %envergadura d'UNA ala = span 
lambda_vec = [0.2];  %taper ratio  
alpha_deg = 4; % Angle d'atac 
U_inf = 27; %[m/s]  (aprox 100km/h)
Nspan = 20; % Nspan discretitza la viga estructural (model 1D en y)
%Malla aerodinamica
Nx = 10;
Ny = 21;
nCases = numel(b_vec)*numel(lambda_vec);

results = repmat(struct('b', [], 'lambda', [], 'S', [], 'AR', [],'L', [], 'CL', [], 'mass', [], 'mu', [], 'CL_over_mu', [], 'sigma_max', [], 'C_sigma', [], 'Csigma_over_mu', [], 'f1', [] ),nCases, 1);
iter = 1;

for i = 1:length(b_vec)
    for j = 1:length(lambda_vec)

        b      = b_vec(i);
        lambda = lambda_vec(j);
        metrics = wing_case(density,youngModulus,poissonRatio, b, lambda, Nspan, Nx, Ny, alpha_deg, U_inf);
        results(iter).y = metrics.y_panel;
        results(iter).b = metrics.b;
        results(iter).lambda = metrics.lambda;
        results(iter).S = metrics.S;
        results(iter).AR = metrics.AR;
        results(iter).L = metrics.L;
        results(iter).L_span = metrics.L_span;
        results(iter).Cp = metrics.Cp;
        results(iter).CL = metrics.CL;
        results(iter).mass = metrics.mass;
        results(iter).mu = metrics.mu;
        results(iter).CL_over_mu = metrics.CL_over_mu;
        results(iter).sigma_max = metrics.sigma_max;
        results(iter).C_sigma = metrics.C_sigma;
        results(iter).Csigma_over_mu = metrics.Csigma_over_mu;
        results(iter).f1 = metrics.f1;

        iter = iter + 1;
    end
end

% figure; 
% hold on; 
% grid on;
% title('Cp distribution for all cases');
% xlabel('y'); ylabel('C_p');
% 
% colors = lines(numel(results));
% 
% for k = 1:numel(results)
%     plot(results(k).y, results(k).Cp, ...
%          'Color', colors(k,:), 'DisplayName', ...
%          sprintf('b=%.1f, λ=%.2f', results(k).b, results(k).lambda));
% end
%legend show;


figure; 
hold on; 
grid on;
title('L distribution for all cases');
xlabel('y'); ylabel('L_{span}');
colors = lines(numel(results));
for k = 1:numel(results)
    plot(results(k).y, results(k).L_span, ...
         'Color', colors(k,:), 'DisplayName', ...
         sprintf('b=%.1f, λ=%.2f', results(k).b, results(k).lambda));
end
legend show;
