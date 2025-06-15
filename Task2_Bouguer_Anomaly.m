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
load([HOME '/Results/gravity_anomaly_mGal.mat'], 'delta_g_mGal')

deltag_mGal = delta_g_mGal;
lmax = 50;                         % Maximum degree/order
R_ref = 2439.4e3;                     % Reference radius (km)
GM = 22031.815e9;              % Mercury GM (km^3/s^2)
G = 6.6743e-11;                     % gravitational constant
rho_crust = 2800;                   % crust density
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
% Using Synthesis
% two_layer_data = load([HOME '/Results/data_bouger_correction_0_50.mat']);
% Model = two_layer_data.Model;
% Data = two_layer_data.data;

%gravity_free_air = deltag_mGal./1e5;
%gravity_free_air(1, 3) = 0;  % Set first row, third column to zero, C_00
%gravity_free_air(4, 3) = 0;  % Set fourth row, third column to zero, C_20
%[free_air_gravity_data] = model_SH_synthesis(lonLim, latLim, height, SHbounds, gravity_free_air, Model1);
%g_freeair = free_air_gravity_data.vec.R;
%g_freeair_mGal = g_freeair * 1e5;  % convert to mGal

% gravity_Model = two_layer_data.V_Model;
% gravity_Model(1,3) = 0;
% gravity_Model(3,3)=0;
% new_gravity_Model = model_SH_synthesis(lonLimT, latLimT, height, SHbounds, gravity_Model, Model);
% bouguer_correction = new_gravity_Model.vec.R;
% g_bouguer_correction_SH_mGal = g_bouguer_correction_SH * 1e5;  % convert to mGal
%BA = g_freeair_mGal - g_bouguer_correction_SH_mGAl;



%Simpler approach similar results
bouguer_correction = 2*pi*G*rho_crust*elevations; % Bouguer correction
bouguer_correction_mGal = bouguer_correction * 1e5; % 1 m/s^2 = 1e5 mGal
scaling_factor = size(bouguer_correction_mGal,1)/size(deltag_mGal,1);
bouguer_correction_mGal = downsize_mean(bouguer_correction_mGal, scaling_factor);
BA = deltag_mGal-bouguer_correction_mGal;

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
set(gca, 'xlim', [0 360]);
set(gca, 'xtick', 0:30:360);
saveas(gcf, 'Figures/Bouguer_Anomaly.svg');

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
set(gca, 'xlim', [0 360]);
set(gca, 'xtick', 0:30:360);
saveas(gcf, 'Figures/Gravity_Anomaly.svg');

figure('units', 'points', 'Position', [0, 0, 455.2441, 0.5*455.2441]);
imagesc(lonT, latT, bouguer_correction_mGal)
c = colorbar;
colormap(broc);
ylabel(c, 'Bouguer Correction (mGal)', 'Interpreter', 'latex', 'Fontsize', aa)
set(gca, 'YDir', 'normal', 'Fontsize', 12)
xlabel('Longitude ($^\circ$)', 'Interpreter', 'latex', 'Fontsize', aa)
ylabel('Latitude ($^\circ$)', 'Interpreter', 'latex', 'Fontsize', aa)
set(gca, 'ylim', [-90 90]);
set(gca, 'ytick', -90:30:90);
set(gca, 'xlim', [0 360]);
set(gca, 'xtick', 0:30:360);
%saveas(gcf, 'Figures/Bouguer_Correction.svg');

figure('units', 'points', 'Position', [0, 0, 455.2441, 0.5*455.2441]);
imagesc(lonT, latT, elevations./1e3)
c = colorbar;
colormap(bam);
ylabel(c, 'Elevation (km)', 'Interpreter', 'latex', 'Fontsize', aa)
set(gca, 'YDir', 'normal', 'Fontsize', 12)
xlabel('Longitude ($^\circ$)', 'Interpreter', 'latex', 'Fontsize', aa)
ylabel('Latitude ($^\circ$)', 'Interpreter', 'latex', 'Fontsize', aa)
set(gca, 'ylim', [-90 90]);
set(gca, 'ytick', -90:30:90);
set(gca, 'xlim', [0 360]);
set(gca, 'xtick', 0:30:360);
%saveas(gcf, 'Figures/Elevation.svg');

% plot_map(lonT, latT, elevations, 'Elevation (m)', 18, 12);