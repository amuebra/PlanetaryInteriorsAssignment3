clear;
close all;
clc;

HOME = pwd;
addpath([HOME '/Data']);
addpath([HOME '/Results']);
addpath([HOME '/Tools'])

% Load the spherical harmonics coefficients from the file
filename = [HOME '/Data/ggmes_50v06_sha.tab'];
maxDegree = 50;
R_ref = 2439.4e3;       % Reference radius in meters
GM = 22031.815e9;       % Mercury GM (m^3/s^2)
coeffs = readmatrix(filename, 'FileType', 'text', 'Delimiter', ',');

% Use coefficients directly as V
V = coeffs;  % Already in [n m Cnm Snm] format
V(3,3) = 0; %set C20 to zero 

% Optional: limit to maxDegree
V = V(V(:,1) <= maxDegree, :);

% Define grid resolution
resolution = 1;
latLimT = [-90+(1/resolution/2) 90-(1/resolution/2) 1/resolution]; 
%lonLimT = [1/resolution/2 360-(1/resolution/2) 1/resolution]     % 0 to 360 degree
lonLimT = [-360/resolution/2  180-(1/resolution/2) 1/resolution]; %-180 to 180 degree

% Generate lat/lon grid
lonT = lonLimT(1):lonLimT(3):lonLimT(2);
%latT = fliplr(latLimT(1):latLimT(3):latLimT(2));
latT = latLimT(1):latLimT(3):latLimT(2);
% LonT = repmat(lonT, length(latT), 1);
% LatT = repmat(latT', 1, length(lonT));

% Define Model structure
Model.Re = R_ref;
Model.GM = GM;
%%
n_vals = unique(V(:,1));
max_n = max(n_vals);
min_n = min(n_vals);

%%
% Spherical harmonic synthesis settings
SHbounds = [1 maxDegree];
height = 0;

% Run synthesis
data = model_SH_synthesis(lonLimT, latLimT, height, SHbounds, V, Model);

% Extract gravity anomaly (in mGal)

delta_g_mGal = data.vec.R * 1e5;

% Plot result
figure
aa = 18;
imagesc(lonT, latT, delta_g_mGal); cc = colorbar;
xlabel('Longitude (\circ)', 'Fontsize', aa)
ylabel('Latitude (\circ)', 'Fontsize', aa)
ylabel(cc, 'Gravity Anomaly (mGal)', 'Fontsize', aa)
set(gca, 'YDir', 'normal', 'Fontsize', aa)

%%
save([HOME '/Results/gravity_anomaly_mGal.mat'], 'delta_g_mGal', 'latT', 'lonT');