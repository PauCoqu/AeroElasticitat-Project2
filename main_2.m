%% Project 2 - Aeroelastic analysis of wing
% Pau Cornudella Quer
% Alex Aleñá 
% 220351 - ADVANCED AEROELASTICITY
%-------------------------------------------
clear; clc; close all

%% STEP 1: Section properties (del pefil) ('beam_properties')  [section]

%Material
density = 2300; %rho Alumini (ENUNCIAT)
youngModulus = 69000000000; % E alumini (ENUNCIAT)
poissonRatio = 0.3; % possion alumini (ENUNCIAT)

%% STEP 2 - Spanwise distributions
%Ala conditions
y0 = 0.34;

% VALORS DE DISSENY FINAL (Projecte 1)
b = 12; %envergadura d'UNA ala = span 
lambda = 0.3;  %taper ratio 

%Section properties (pel beam_section_properties)
geometry = "naca"; %Perfil Naca (ENUNCIAT) "naca"

%param(1)= chord ; param(2) = thickness (0.1*chord per al NACA 0010)
%param(3) = skin thickness
%(x1,t1) = (lenght,thickness) = (0.22,0.02)   ->  spar 1
%(x2,t2) = (lenght,thickness) = (0.52,0.02)   ->  spar 2

         % p1   p2    p3   x1   t1   x2   t2 
param = [0.85,0.085];%,0.015,0.22,0.02,0.52,0.02]; 

alpha_deg = 4; % Angle d'atac 
U_inf = 27 ;%[27 110 130]; %[m/s]  (aprox 100km/h)

% Malla estructural i Malla aerodinamica
Nx = 20; %Numero d'elements
Ny = 20;

%Modal reduction
i_modes = [1:4];

metrics = wing_case(geometry, param, density,youngModulus,poissonRatio,b,lambda,Nx, Ny,U_inf,y0,i_modes);
