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

alpha_deg = 4; % Angle d'atac 
U_inf = 27 ;%[27 110 130]; %[m/s]  (aprox 100km/h)

% Malla estructural i Malla aerodinamica
Nx = 20; %Numero d'elements
Ny = 20;

metrics = wing_case(density,youngModulus,poissonRatio,b,lambda,Nx, Ny,U_inf,y0);
