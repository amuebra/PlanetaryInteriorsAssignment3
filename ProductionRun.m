% main file for the complete GSH circle for a particular model
clear;
close all;
clc;

HOME = pwd;
addpath([HOME '/Data']);
addpath([HOME '/Results'])
addpath([HOME '/Tools']);

% Model Construction

new_model = 1;

if new_model == 1

  % Construct new model
  
  Model = struct();
  
  Model.number_of_layers = 2;
  Model.name = 'two_layer_planet';
  
  % Additional variables
  Model.GM = 22031.815E9; % in m^3/s^2
  Model.Re = 2439.4E3; %Reference radius in m
  Model.geoid = 'none';
  Model.nmax = 179;     
  Model.correct_depth = 0;
  
  % Top layer
  %Topography (maxtrix latitude and logitude in m) 
  %Model.l1.bound = gmt2matrix(load([HOME '/Data/crust1.bd1.gmt'])).*1e3;  % meters with respect to reference sphere
  tmp = load([HOME '/Results/Airy_thickness.mat'], 'crust_thickness','longitudes', 'latitudes');
  crust_thickness = tmp.crust_thickness;
  longitudes = tmp.longitudes;
  latitudes = tmp.latitudes;
  Model.l1.bound = crust_thickness;
  Model.l1.dens  = 2800;
  
  %% Second layer
  %Model.l2.bound = -50000+Model.l1.bound;     % meters with respect to reference sphere ( or -50000.*ones(size(TOPO)) make a loop around Second layer fr M1 and M2
  Model.l2.bound = -50000.*ones(size(crust_thickness));
  Model.l2.dens  = 3200;	   % Density in kg/m3
  
  % Bottom bound
  Model.l3.bound = -100000;    % meters with respect to reference sphere (make it large enough, should be lower
  
  % Save model in .mat file for use of the new software
  
  save([HOME '/Data/' Model.name '.mat'],'Model')

else
  % Load previous saved model

  model_name = 'two_layer_planet';
  load([HOME '/Data/' Model.name '.mat']);
end

%%%%%%%%%%%%%%%%%%% Computation area %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Part that can be modified %%%%%%%%%%%%%%%%%%%%%%%

latLim =    [-89.5 89.5 1];  % [deg] min latitude, max latitude, resolution latitude (preferable similar to latitude)
lonLim =    [-180 180 1];% [deg] min longitude, max longitude, resolution longitude (preferable similar to latitude)
height =    0; % height of computation above spheroid
SHbounds =  [0 179]; % Truncation settings: lower limit, upper limit SH-coefficients used

%%%%%%%%%%%%%% Part that can be modified %%%%%%%%%%%%%%%%%%%%%%%

%% Global Spherical Harmonic Analysis 

tic;
[V_Model] = segment_2layer_model(Model.l1.bound,Model.l2.bound,Model.l3.bound,Model.l1.dens,Model.l2.dens,62000,Model); % order, number, S-coefficentes, c-coefficents
%or something similar it shows spherical harmonics
toc

%% Global Spherical Harmonic Synthesis

tic;
[data] = model_SH_synthesis(lonLim,latLim,height,SHbounds,V_Model,Model);
toc

%% Save data

save([HOME '/Results/data_' Model.name '_' num2str(SHbounds(1)) '_' num2str(SHbounds(2)) '.mat'],'data','V_Model','Model')
% Extract gravity anomaly (in mGal)
delta_g_mGal = data.vec.R * 1e5;

%% plot results
figure
aa = 18;
imagesc(lonLim, latLim, delta_g_mGal);
c = colorbar;
ylabel(c, 'Gravity Anomaly (mGal)', 'Interpreter', 'latex', 'Fontsize', aa)
set(gca, 'YDir', 'normal', 'Fontsize', 12)
xlabel('Longitude ($^\circ$)', 'Interpreter', 'latex', 'Fontsize', aa)
ylabel('Latitude ($^\circ$)', 'Interpreter', 'latex', 'Fontsize', aa)
set(gca, 'ylim', [-90 90]);
set(gca, 'ytick', -90:30:90);
set(gca, 'xlim', [-180 180]);
set(gca, 'xtick', -180:30:180);
