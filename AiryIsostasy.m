% Mercury_Isostasy.m
% Airy isostasy model for Mercury's crust-mantle interface
clear;
close all;
clc;

HOME = pwd;
addpath([HOME '/Data']);
addpath([HOME '/Results']);

% --- Load Topography Data ---
load([HOME '/Results/elevations.mat'], 'elevations', 'longitudes', 'latitudes')  % Expects `elevations` defined (in meters)

% %% --- Downsample topography early to save memory ---
% h = elevations;
% h = double(h);
% step = 20;  % Keep every 20th point (adjust to fit your memory)
% h = h(1:step:end, 1:step:end);  % Keep every 100th row and column

%% --- Physical Parameters (adjust as needed) ---
D = 35e3;            % Crustal thickness in meters (e.g., 35 km)
rho_crust = 2800;     % Crustal density (kg/m^3)
rho_mantle = 3200;    % Mantle density (kg/m^3)
h = elevations;

%% --- Compute compensation at Moho using Airy isostasy ---
% r_continent is the variation at the crust-mantle interface (Moho)
r_continent = h .* (rho_crust / (rho_mantle - rho_crust));
% --- Compute actual crustal thickness ---
crust_thickness = D + r_continent;

%% --- Plot the results ---
figure;
imagesc(longitudes, latitudes, crust_thickness);  % Automatically scales axes
set(gca, 'YDir', 'normal');          % Flip Y so North is up
xlabel('Longitude (°E)');
ylabel('Latitude (°)');
title('Mercury Crustal Thickness from Airy Isostasy');
colorbar;
colormap(jet);
axis equal tight;

% --- Optional: Save results ---
save([HOME '/Results/Airy_r_contintent.mat'], 'r_continent', 'longitudes', 'latitudes');
save([HOME '/Results/Airy_thickness.mat'], 'crust_thickness', 'longitudes', 'latitudes');
