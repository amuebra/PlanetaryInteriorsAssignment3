%%
clear
close all
clc

HOME = pwd;
addpath([HOME '/Data']);
addpath([HOME '/Results']);
addpath([HOME '/Tools']);

%% parameters
rho_crust = 2850;
rho_mantle = 3100;
R_crust = 2439.36;
R_mantle = 2410;
D = (R_crust - R_mantle) * 1e3;

%% load data
% load('Data/acton.mat')
data = read(Tiff('Mercury_Messenger_USGS_DEM_Global_665m_v2.tif'));
%data2 = read(Tiff('MercuryMessengerUSGS_MAP2_EQUI.tif'));

%% processing
resolution = 64;
latitudes = -90 : 1/resolution : 90;
latitudes = latitudes(1:end-1);
longitudes = -180 : 1/resolution : 180;
longitudes = longitudes(1:end-1);

multiplier = 0.5;
elevations = multiplier * double(data);
h = elevations;

%% calculate root of topography
r = (h*rho_crust)/(rho_mantle - rho_crust); % changes of the mantle crust boundary (root)
t = D - r;

%% plotting
figure;
imagesc(longitudes, latitudes, t)
colorbar
% colormap(acton)

%% analysis
disp(D)

% save('r.mat', 'r')
% save('elevations.mat', 'elevations');