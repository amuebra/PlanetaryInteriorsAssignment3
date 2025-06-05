%%
clear
close all
clc

%% constants
G = 6.67430e-11;
rho_crust = 2800;

%% load data
% load('Data/acton.mat')
data = read(Tiff('Mercury_Messenger_USGS_DEM_Global_665m_v2.tif'));

%% processing
resolution = 64;
latitudes = -90 : 1/resolution : 90;
latitudes = latitudes(1:end-1);
longitudes = -180 : 1/resolution : 180;
longitudes = longitudes(1:end-1);

multiplier = 0.5;
elevations = multiplier * double(data);

%% calculate Bouguer correction
deltag_b = 2*pi*G*rho_crust*elevations; % Bouguer correction

%% plotting
% topography = elevations/1000;
% imagesc(longitudes, latitudes, topography)
% colorbar
% colormap(acton)

imagesc(longitudes, latitudes, deltag_b)
colorbar

%% analysis
disp(max(max(data./2)))
disp(min(min(data./2)))
