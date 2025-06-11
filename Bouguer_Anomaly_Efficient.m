clear;
close all;
clc;

HOME = pwd;
addpath([HOME '/Data']);
addpath([HOME '/Results']);
addpath([HOME '/Tools']);
addpath([HOME '/ScientificColourMaps8/vik']);
addpath([HOME '/ScientificColourMaps8/cork']);
addpath([HOME '/ScientificColourMaps8/broc']);
addpath([HOME '/ScientificColourMaps8/bam']);
load('vik.mat');
load('cork.mat');
load('broc.mat');
load('bam.mat');

%% Parameters
load([HOME '/Results/elevations.mat'], 'elevations')
load([HOME '/Results/gravity_anomaly_mGal.mat'], 'delta_g_mGal', 'latT', 'lonT')
deltag_mGal = delta_g_mGal;
lmax = 100;                         % Maximum degree/order
R_ref = 2439.4;                     % Reference radius (km)
GM = 22031.8150000000;              % Mercury GM (km^3/s^2)
G = 6.6743e-11;                     % gravitational constant
rho_crust = 2800;                   % crust density
resolution = 4;                     % degrees (1 = 1x1°, 4 = 0.25x0.25°) from GravityChangedLatitude.m

%% CALCULATE BOUGUER ANOMALY
% -------------------------------------------------------------------------
deltag_b = 2*pi*G*rho_crust*elevations; % Bouguer correction
deltag_b_mGal = deltag_b * 1e5; % 1 m/s^2 = 1e5 mGal
scaling_factor = size(deltag_b_mGal,1)/size(deltag_mGal,1);

deltag_b_mGal = downsize_mean(deltag_b_mGal, scaling_factor);
BA = deltag_mGal - deltag_b_mGal; % Bouguer Anomaly in mGal

% -------------------------------------------------------------------------
%% PLOTTING
% -------------------------------------------------------------------------
aa = 18;
figure('units', 'points', 'Position', [0, 0, 455.2441, 0.5*455.2441]);
imagesc(lonT, latT, BA)
c = colorbar;
colormap(vik);
ylabel(c, 'Bouguer Anomaly (mGal)', 'Interpreter', 'latex', 'Fontsize', aa)
set(gca, 'YDir', 'normal', 'Fontsize', 12)
xlabel('Longitude ($^\circ$)', 'Interpreter', 'latex', 'Fontsize', aa)
ylabel('Latitude ($^\circ$)', 'Interpreter', 'latex', 'Fontsize', aa)
set(gca, 'ylim', [-90 90]);
set(gca, 'ytick', -90:30:90);
set(gca, 'xlim', [-180 180]);
set(gca, 'xtick', -180:30:180);
saveas(gcf, 'figures/Bouguer_Anomaly.svg');

figure('units', 'points', 'Position', [0, 0, 455.2441, 0.5*455.2441]);
imagesc(lonT, latT, deltag_mGal)
c = colorbar;
colormap(cork);
ylabel(c, 'Free-air Gravity Anomaly (mGal)', 'Interpreter', 'latex', 'Fontsize', aa)
set(gca, 'YDir', 'normal', 'Fontsize', 12)
xlabel('Longitude ($^\circ$)', 'Interpreter', 'latex', 'Fontsize', aa)
ylabel('Latitude ($^\circ$)', 'Interpreter', 'latex', 'Fontsize', aa)
set(gca, 'ylim', [-90 90]);
set(gca, 'ytick', -90:30:90);
set(gca, 'xlim', [-180 180]);
set(gca, 'xtick', -180:30:180);
saveas(gcf, 'figures/Gravity_Anomaly.svg');

figure('units', 'points', 'Position', [0, 0, 455.2441, 0.5*455.2441]);
imagesc(lonT, latT, deltag_b_mGal)
c = colorbar;
colormap(broc);
ylabel(c, 'Bouguer Correction (mGal)', 'Interpreter', 'latex', 'Fontsize', aa)
set(gca, 'YDir', 'normal', 'Fontsize', 12)
xlabel('Longitude ($^\circ$)', 'Interpreter', 'latex', 'Fontsize', aa)
ylabel('Latitude ($^\circ$)', 'Interpreter', 'latex', 'Fontsize', aa)
set(gca, 'ylim', [-90 90]);
set(gca, 'ytick', -90:30:90);
set(gca, 'xlim', [-180 180]);
set(gca, 'xtick', -180:30:180);
saveas(gcf, 'figures/Bouguer_Correction.svg');

figure('units', 'points', 'Position', [0, 0, 455.2441, 0.5*455.2441]);
imagesc(lonT, latT, elevations)
c = colorbar;
colormap(bam);
ylabel(c, 'Elevation (m)', 'Interpreter', 'latex', 'Fontsize', aa)
set(gca, 'YDir', 'normal', 'Fontsize', 12)
xlabel('Longitude ($^\circ$)', 'Interpreter', 'latex', 'Fontsize', aa)
ylabel('Latitude ($^\circ$)', 'Interpreter', 'latex', 'Fontsize', aa)
set(gca, 'ylim', [-90 90]);
set(gca, 'ytick', -90:30:90);
set(gca, 'xlim', [-180 180]);
set(gca, 'xtick', -180:30:180);
saveas(gcf, 'figures/Elevation.svg');

% plot_map(lonT, latT, elevations, 'Elevation (m)', 18, 12);