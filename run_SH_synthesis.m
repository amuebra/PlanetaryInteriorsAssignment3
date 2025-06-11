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
coeffs = [[0,0,0,0,0,0]; coeffs];

% Use coefficients directly as V
V = coeffs;  % Already in [n m Cnm Snm] format
V(4,3) = 0; %set C20 to zero 

% Optional: limit to maxDegree
V = V(V(:,1) <= maxDegree, :);

% Define grid resolution
resolution = 1;
latLimT = [-90+(1/resolution/2) 90-(1/resolution/2) 1/resolution]; 
%lonLimT = [1/resolution/2 360-(1/resolution/2) 1/resolution]     % 0 to 360 degree
lonLimT = [-180+(1/resolution/2) 180-(1/resolution/2) 1/resolution]; %-180 to 180 degree

% Generate lat/lon grid
lonT = lonLimT(1):lonLimT(3):lonLimT(2);
%latT = fliplr(latLimT(1):latLimT(3):latLimT(2));
latT = latLimT(1):latLimT(3):latLimT(2);

%% perform computation
% Define Model structure
Model.Re = R_ref;
Model.GM = GM;

% Spherical harmonic synthesis settings
SHbounds = [1 maxDegree];
height = 0;

% Run synthesis
data = model_SH_synthesis(lonLimT, latLimT, height, SHbounds, V, Model);

% Extract gravity anomaly (in mGal)
delta_g_mGal = data.vec.R * 1e5;


%% plot results
figure
aa = 18;
imagesc(lonT, latT, delta_g_mGal);
c = colorbar;
colormap(turbo)
ylabel(c, 'Gravity Anomaly (mGal)', 'Interpreter', 'latex', 'Fontsize', aa)
set(gca, 'YDir', 'normal', 'Fontsize', 12)
xlabel('Longitude ($^\circ$)', 'Interpreter', 'latex', 'Fontsize', aa)
ylabel('Latitude ($^\circ$)', 'Interpreter', 'latex', 'Fontsize', aa)
clim(-80, 120),
% set(gca, 'ylim', [-90 90]);
% set(gca, 'ytick', -90:30:90);
% set(gca, 'xlim', [-180 180]);
% set(gca, 'xtick', -180:30:180);

%% save data
save([HOME '/Results/gravity_anomaly_mGal.mat'], 'delta_g_mGal', 'latT', 'lonT');