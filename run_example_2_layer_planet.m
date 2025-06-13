% main file for the complete GSH circle for a particular model
clear;
close all;
clc;

HOME = pwd;
addpath([HOME '/Data']);
addpath([HOME '/Tools']);
addpath([HOME '/Results']);
load([HOME '/Results/elevations.mat'], 'elevations')
% Model Construction

new_model = 1;

if new_model == 1

  % Construct new model
  
  Model = struct();
  
  Model.number_of_layers = 2;
  Model.name = 'bouger_correction';
  
  % Additional variables
  Model.GM = 22031.815E9; % in m^3/s^2
  Model.Re = 2439.4E3; %Reference radius in m
  Model.geoid = 'none';
  Model.nmax = 50;     
  Model.correct_depth = 0;
  
  % Top layer
  %Topography (maxtrix latitude and logitude in m) 
  Model.l1.bound = elevations;  % meters with respect to reference sphere 
  Model.l1.dens  = 2800;
  
  % Second layer
  Model.l2.bound = -50000*ones(size(elevations));     % meters with respect to reference sphere ( or -50000.*ones(size(TOPO)) make a loop around Second layer fr M1 and M2
  Model.l2.dens  = 3200;	   % Density in kg/m3
  
  % Bottom bound
  Model.l3.bound = -100000;    % meters with respect to reference sphere (make it large enough, should be lower
  
  % Save model in .mat file for use of the new software
  
  save([HOME '/Data/' Model.name '.mat'],'Model')

else
  % Load previous saved model

  model_name = 'bouger_correction';
  load([HOME '/Data/' Model.name '.mat']);
end

%%%%%%%%%%%%%%%%%%% Computation area %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Part that can be modified %%%%%%%%%%%%%%%%%%%%%%%

latLim =    [-89.5 89.5 1];  % [deg] min latitude, max latitude, resolution latitude (preferable similar to latitude)
lonLim =    [0.05 359.5 1];% [deg] min longitude, max longitude, resolution longitude (preferable similar to latitude)
height =    0; % height of computation above spheroid
SHbounds =  [0 50]; % Truncation settings: lower limit, upper limit SH-coefficients used

%%%%%%%%%%%%%% Part that can be modified %%%%%%%%%%%%%%%%%%%%%%%

%% Global Spherical Harmonic Analysis 

tic;
[V_Model] = segment_2layer_model(Model.l1.bound,Model.l2.bound,Model.l3.bound,Model.l1.dens,Model.l2.dens,25000,Model); % order, number, S-coefficentes, c-coefficents
%or something similar it shows spherical harmonics
toc

%% Global Spherical Harmonic Synthesis

tic;
[data] = model_SH_synthesis(lonLim,latLim,height,SHbounds,V_Model,Model);
toc

%% Save data

save([HOME '/Results/data_' Model.name '_' num2str(SHbounds(1)) '_' num2str(SHbounds(2)) '.mat'],'data','V_Model','Model
