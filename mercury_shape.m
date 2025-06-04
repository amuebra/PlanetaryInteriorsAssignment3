%%
clear
close all
clc

%% load data
% load('Data/acton.mat')
data = read(Tiff('Mercury_Messenger_USGS_DEM_Global_665m_v2.tif'));
data2 = read(Tiff('MercuryMessengerUSGS_MAP2_EQUI.tif'));

%% processing
resolution = 64;
latitudes = -90 : 1/resolution : 90;
latitudes = latitudes(1:end-1);
longitudes = -180 : 1/resolution : 180;
longitudes = longitudes(1:end-1);

multiplier = 0.5;
elevations = multiplier * data;

%% plotting
imagesc(longitudes, latitudes, elevations)
colorbar
% colormap(acton)

%% analysis
disp(max(max(data./2)))
disp(min(min(data./2)))
