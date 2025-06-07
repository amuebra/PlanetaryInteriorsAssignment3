% Mercury_Isostasy.m
% Airy isostasy model for Mercury's crust-mantle interface
clear;
close all;
clc;

HOME = pwd;
addpath([HOME '/Data']);
addpath([HOME '/Results']);

% --- Load Topography Data ---
load([HOME '/Results/elevations.mat'], 'elevations')  % Expects `elevations` defined (in meters)
h = elevations;

% --- Physical Parameters (adjust as needed) ---
D = 35000;            % Crustal thickness in meters (e.g., 35 km)
rho_crust = 2800;     % Crustal density (kg/m^3)
rho_mantle = 3300;    % Mantle density (kg/m^3)

% --- Compute compensation at Moho using Airy isostasy ---
% r_continent is the variation at the crust-mantle interface (Moho)
r_continent = h .* (rho_crust / (rho_mantle - rho_crust)) - D;

% Optional: plot the results
figure;
subplot(2,1,1);
plot(h, 'b');
title('Surface Topography');
ylabel('Elevation (m)');

subplot(2,1,2);
plot(r_continent, 'r');
title('Crust-Mantle Interface (Airy Isostasy)');
ylabel('Moho Depth Variation (m)');
xlabel('Sample Index');

% --- Optional: Save results ---
%save('Moho_Depth_Mercury.mat', 'r_continent');
