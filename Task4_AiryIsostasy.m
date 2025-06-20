% Mercury_Isostasy.m
% Airy isostasy model for Mercury's crust-mantle interface
clear;
close all;
clc;

HOME = pwd;
addpath([HOME '/Data']);
addpath([HOME '/Results']);
load('vik.mat');

% --- Load Topography Data ---
load([HOME '/Results/elevations.mat'], 'elevations', 'longitudes', 'latitudes')  % Expects `elevations` defined (in meters)


%% --- Physical Parameters (adjust as needed) ---
D = 36459;            % Crustal thickness in meters (e.g., 35 km)
rho_crust = 2800;     % Crustal density (kg/m^3)
rho_mantle = 3200;    % Mantle density (kg/m^3)
h = elevations;

%% --- Compute compensation at Moho using Airy isostasy ---
% r_continent is the variation at the crust-mantle interface (Moho)
r_continent = h .* (rho_crust / (rho_mantle - rho_crust));
% --- Compute actual crustal thickness ---
crust_thickness = D + r_continent;

%% --- Plot the results ---
aa =16;
figure;
imagesc(longitudes, latitudes, crust_thickness./1e3);  % Automatically scales axes
set(gca, 'YDir', 'normal','Fontsize', aa);          % Flip Y so North is up
xlabel('Longitude ($^\circ$)', 'Interpreter', 'latex', 'Fontsize', aa)
ylabel('Latitude ($^\circ$)', 'Interpreter', 'latex', 'Fontsize', aa)
%title('Mercury Crustal Thickness from Airy Isostasy');
c = colorbar;
ylabel(c, 'Airy Crust Thickness (km)', 'Interpreter', 'latex', 'Fontsize', aa);
c = colormap(vik);
set(gca, 'ylim', [-90 90]);
set(gca, 'ytick', -90:30:90);
set(gca, 'xlim', [0 360]);
set(gca, 'xtick', 0:30:360);
axis equal tight;
saveas(gcf, 'Figures/Airy_Thickness.svg');
% --- Optional: Save results ---
%save([HOME '/Results/Airy_r_contintent.mat'], 'r_continent', 'longitudes', 'latitudes');
save([HOME '/Results/Airy_thickness.mat'], 'crust_thickness', 'longitudes', 'latitudes');
