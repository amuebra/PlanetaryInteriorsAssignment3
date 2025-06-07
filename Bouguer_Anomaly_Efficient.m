clear;
close all;
clc;

HOME = pwd;
addpath([HOME '/Data']);
addpath([HOME '/Results'])

%% Parameters
load([HOME '/Results/elevations.mat'], 'elevations')
load([HOME '/Results/r4_delta_g_mGal.mat'], 'delta_g_mGal')
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
[a, b] = size(deltag_b_mGal);
if mod(a, scaling_factor) ~= 0 || mod(b, scaling_factor) ~= 0
    error('Matrix dimensions must be divisible by the reduction factor.');
end

% Reshape and average to reduce size of deltag_b
deltag_b_mGal = mean(reshape(deltag_b_mGal, scaling_factor, a/scaling_factor, scaling_factor, b/scaling_factor), [1 3]);
deltag_b_mGal = squeeze(deltag_b_mGal);
BA = deltag_mGal - deltag_b_mGal; % Bouguer Anomaly in mGal

figure;
imagesc(longitudes, latitudes, BA)
c = colorbar;
c.Label.String = 'Bouguer Anomaly in mGal';