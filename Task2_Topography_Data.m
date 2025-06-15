%%
clear
close all
clc

HOME = pwd;
addpath([HOME '/Data']);
addpath([HOME '/Results']);
addpath([HOME '/Tools']);

%% load data
% load('Data/acton.mat')
data = read(Tiff([HOME '/Data/Mercury_Messenger_USGS_DEM_Global_665m_v2.tif']));
%data2 = read(Tiff('MercuryMessengerUSGS_MAP2_EQUI.tif'));

%% processing
resolution = 64;
latitudes = -90 : 1/resolution : 90;
latitudes = latitudes(1:end-1);
longitudes = 0 : 1/resolution : 360;
longitudes = longitudes(1:end-1);

multiplier = 0.5;
elevations = multiplier * double(data);
elevations = flipud(elevations);
%elevations = Europe_centered(elevations);

%% plotting
imagesc(longitudes, latitudes, elevations);
colorbar;
set(gca,'YDir','normal');
set(gca, 'ylim', [-90 90]);
set(gca, 'ytick', -90:30:90);
% colormap(acton)

%% analysis
%disp(max(max(data./2)))
%disp(min(min(data./2)))

elevations = downsize_mean(elevations, 64);
longitudes = downsize_mean(longitudes, 64);
latitudes = downsize_mean(latitudes, 64);
save([HOME '/Results/elevations.mat'], 'elevations', 'longitudes', 'latitudes');