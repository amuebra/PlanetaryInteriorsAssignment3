clear;
close all;
clc;

HOME = pwd;
addpath([HOME '/Data']);
addpath([HOME '/Results'])
addpath([HOME '/Tools']);

%% Parameters
load([HOME '/Results/elevations.mat'], 'elevations')
load([HOME '/Results/gravity_anomaly_mGal.mat'], 'delta_g_mGal', 'lat', 'lon')
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
figure;
imagesc(lon, lat, BA)
c = colorbar;
c.Label.String = 'Bouguer Anomaly in mGal';
set(gca,'YDir','normal')
set(gca, 'ylim', [-90 90]);
set(gca, 'ytick', -90:30:90);

figure;
imagesc(lon, lat, deltag_mGal)
c = colorbar;
c.Label.String = 'Gravity Anomaly in mGal';
set(gca,'YDir','normal')
set(gca, 'ylim', [-90 90]);
set(gca, 'ytick', -90:30:90);

figure;
imagesc(lon, lat, deltag_b_mGal)
c = colorbar;
c.Label.String = 'Bouguer Correction in mGal';
set(gca,'YDir','normal')
set(gca, 'ylim', [-90 90]);
set(gca, 'ytick', -90:30:90);

figure;
imagesc(lon, lat, elevations)
c = colorbar;
c.Label.String = 'elevation in m';
set(gca,'YDir','normal')
set(gca, 'ylim', [-90 90]);
set(gca, 'ytick', -90:30:90);