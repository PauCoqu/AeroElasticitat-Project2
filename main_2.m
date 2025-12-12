%% Project 2 - Aeroelastic analysis of wing
% Pau Cornudella Quer
% Alex Aleñá 
% 220351 - ADVANCED AEROELASTICITY
%-------------------------------------------
clear; clc; close all

%% PROJECTE 1:

%% STEP 1: Section properties (del pefil) ('beam_properties')  [section]

%Material
density = 2300; %rho Alumini (ENUNCIAT)
youngModulus = 69000000000; % E alumini (ENUNCIAT)
poissonRatio = 0.3; % possion alumini (ENUNCIAT)

%% STEP 2 - Spanwise distributions

% VALORS DE DISSENY FINAL (Projecte 1)
b = 12 ; %envergadura d'UNA ala = span 
lambda = 0.2;  %taper ratio  


alpha_deg = 4; % Angle d'atac 
U_inf = 27; %[m/s]  (aprox 100km/h)
Nspan = 20; % Nspan discretitza la viga estructural (model 1D en y)
%Malla aerodinamica
Nx = 10;
Ny = 21;

metrics = wing_case(density,youngModulus,poissonRatio, b, lambda, Nspan, Nx, Ny, alpha_deg, U_inf);

%% 

x_panel = metrics.x_panel;
y_panel = metrics.y_panel;
c_root  = metrics.c_root;

M = 0;   % incompressible
k = 0;   % quasi-steady

AIC0 = buildAIC(x_panel, y_panel, c_root, M, k);

fprintf("AIC0 size = %dx%d\n", size(AIC0,1), size(AIC0,2));

