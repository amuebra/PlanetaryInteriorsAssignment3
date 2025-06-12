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
%load('vik.mat');
%load('cork.mat');
%load('broc.mat');
%load('bam.mat');

%% Parameters
load([HOME '/Results/elevations.mat'], 'elevations')
load([HOME '/Results/gravity_anomaly_mGal.mat'], 'delta_g_mGal')

deltag_mGal = delta_g_mGal;
lmax = 100;                         % Maximum degree/order
R_ref = 2439.4;                     % Reference radius (km)
GM = 22031.8150000000;              % Mercury GM (km^3/s^2)
G = 6.6743e-11;                     % gravitational constant
rho_crust = 2500;                   % crust density
resolution = 1;                     % degrees (1 = 1x1°, 4 = 0.25x0.25°) from GravityChangedLatitude.m
height = 0;
SHbounds = [1 50];

%
latLimT = [-90+(1/resolution/2) 90-(1/resolution/2) 1/resolution]; 
lonLimT = [1/resolution/2 360-(1/resolution/2) 1/resolution];
latT = latLimT(1):latLimT(3):latLimT(2);
lonT = lonLimT(1):lonLimT(3):lonLimT(2);

%% CALCULATE BOUGUER ANOMALY
% -------------------------------------------------------------------------
two_layer_data = load([HOME '/Results/data_bouger_correction_0_50.mat']);
Model = two_layer_data.Model;
gravity_Model = two_layer_data.V_Model;
gravity_Model(1,3) = 0;
gravity_Model(3,3)=0;
new_gravity_Model = model_SH_synthesis(lonLimT, latLimT, height, SHbounds, gravity_Model, Model);
bouger_correction = new_gravity_Model.vec.R;
bouger_correction_mGal = bouger_correction *1e5;
%deltag_b = 2*pi*G*rho_crust*elevations; % Bouguer correction
%deltag_b_correction = deltag_b-deltag_correction;

%deltag_b_mGal = deltag_b * 1e5; % 1 m/s^2 = 1e5 mGal
%scaling_factor = size(deltag_b_mGal,1)/size(deltag_mGal,1);

%deltag_b_mGal = downsize_mean(deltag_b_mGal, scaling_factor);
BA = deltag_mGal - bouger_correction_mGal; % Bouguer Anomaly in mGal

% -------------------------------------------------------------------------
%% PLOTTING
% -------------------------------------------------------------------------
aa = 18;
figure('units', 'points', 'Position', [0, 0, 455.2441, 0.5*455.2441]);
imagesc(lonT, latT, BA)
c = colorbar;
%colormap(vik);
ylabel(c, 'Bouguer Anomaly (mGal)', 'Interpreter', 'latex', 'Fontsize', aa)
set(gca, 'YDir', 'normal', 'Fontsize', 12)
xlabel('Longitude ($^\circ$)', 'Interpreter', 'latex', 'Fontsize', aa)
ylabel('Latitude ($^\circ$)', 'Interpreter', 'latex', 'Fontsize', aa)
set(gca, 'ylim', [-90 90]);
set(gca, 'ytick', -90:30:90);
%set(gca, 'xlim', [-180 180]);
%set(gca, 'xtick', -180:30:180);
saveas(gcf, 'figures/Bouguer_Anomaly.svg');

figure('units', 'points', 'Position', [0, 0, 455.2441, 0.5*455.2441]);
imagesc(lonT, latT, deltag_mGal)
c = colorbar;
%colormap(cork);
ylabel(c, 'Free-air Gravity Anomaly (mGal)', 'Interpreter', 'latex', 'Fontsize', aa)
set(gca, 'YDir', 'normal', 'Fontsize', 12)
xlabel('Longitude ($^\circ$)', 'Interpreter', 'latex', 'Fontsize', aa)
ylabel('Latitude ($^\circ$)', 'Interpreter', 'latex', 'Fontsize', aa)
set(gca, 'ylim', [-90 90]);
set(gca, 'ytick', -90:30:90);
%set(gca, 'xlim', [-180 180]);
%set(gca, 'xtick', -180:30:180);
saveas(gcf, 'figures/Gravity_Anomaly.svg');

figure('units', 'points', 'Position', [0, 0, 455.2441, 0.5*455.2441]);
imagesc(lonT, latT, bouger_correction_mGal)
c = colorbar;
%colormap(broc);
colormap(turbo);
ylabel(c, 'Bouguer Correction (mGal)', 'Interpreter', 'latex', 'Fontsize', aa)
set(gca, 'YDir', 'normal', 'Fontsize', 12)
xlabel('Longitude ($^\circ$)', 'Interpreter', 'latex', 'Fontsize', aa)
ylabel('Latitude ($^\circ$)', 'Interpreter', 'latex', 'Fontsize', aa)
set(gca, 'ylim', [-90 90]);
set(gca, 'ytick', -90:30:90);
%set(gca, 'xlim', [-180 180]);
%set(gca, 'xtick', -180:30:180);
%saveas(gcf, 'figures/Bouguer_Correction.svg');

figure('units', 'points', 'Position', [0, 0, 455.2441, 0.5*455.2441]);
imagesc(lonT, latT, elevations)
c = colorbar;
colormap(turbo);
ylabel(c, 'Elevation (m)', 'Interpreter', 'latex', 'Fontsize', aa)
set(gca, 'YDir', 'normal', 'Fontsize', 12)
xlabel('Longitude ($^\circ$)', 'Interpreter', 'latex', 'Fontsize', aa)
ylabel('Latitude ($^\circ$)', 'Interpreter', 'latex', 'Fontsize', aa)
set(gca, 'ylim', [-90 90]);
set(gca, 'ytick', -90:30:90);
%set(gca, 'xlim', [-180 180]);
%set(gca, 'xtick', -180:30:180);
%saveas(gcf, 'figures/Elevation.svg');

% plot_map(lonT, latT, elevations, 'Elevation (m)', 18, 12);