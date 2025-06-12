% main file for the complete GSH circle for a particular model
clear;
close all;
clc;


HOME = pwd;
addpath([HOME '/Data']);
addpath([HOME '/Results'])
addpath([HOME '/Tools']);

% Load Data
tpl = load([HOME '/Results/elevations.mat'], 'elevations');
elevations = tpl.elevations;
amp = load([HOME '/Results/Airy_thickness.mat'], 'crust_thickness','longitudes', 'latitudes');
crust_thickness = amp.crust_thickness;
longitudes = amp.longitudes;
latitudes = amp.latitudes;
fmp = load([HOME '/Results/Flexural_thickness.mat'], 'mapf');
flexural_crust_thickness = fmp.mapf;

% Model Construction

new_model = 1;

if new_model == 1

  % Construct new model
  
  Model = struct();
  
  Model.number_of_layers = 2;
  Model.name = 'Airy_variance';
  
  % Additional variables
  Model.GM = 22031.815E9; % in m^3/s^2
  Model.Re = 2439.4E3; %Reference radius in m
  Model.geoid = 'none';
  Model.nmax = 179;     
  Model.correct_depth = 0;
  
  % Top layer
  %Topography (maxtrix latitude and logitude in m) 
  %Model.l1.bound = gmt2matrix(load([HOME '/Data/crust1.bd1.gmt'])).*1e3;  % meters with respect to reference sphere
  Model.l1.bound = elevations;
  Model.l1.dens  = 2800;
  
  %% Second layer
  %Model.l2.bound = -50000+Model.l1.bound;     % meters with respect to reference sphere ( or -50000.*ones(size(TOPO)) make a loop around Second layer fr M1 and M2
  Model.l2.bound = -crust_thickness;
  Model.l2.dens  = 3200;	   % Density in kg/m3
  
  % Bottom bound
  Model.l3.bound = -100000;    % meters with respect to reference sphere (make it large enough, should be lower
  
  % Save model in .mat file for use of the new software
  
  save([HOME '/Data/' Model.name '.mat'],'Model')

else
  % Load previous saved model

  model_name = 'Airy_variance';
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

%% Second model with modified layer 2

Model_B = Model;  % Copy the base model
Model_B.name = 'Flexural_variance';

% Modify the second layer boundary (example: 10 km deeper)
Model_B.l2.bound = -flexural_crust_thickness;

% Run SH analysis and synthesis
[V_Model_B] = segment_2layer_model(Model_B.l1.bound, Model_B.l2.bound, Model_B.l3.bound, ...
                                   Model_B.l1.dens, Model_B.l2.dens, 62000, Model_B);

[data_B] = model_SH_synthesis(lonLim, latLim, height, SHbounds, V_Model_B, Model_B);

% Save second model and results
save([HOME '/Results/data_' Model_B.name '_' num2str(SHbounds(1)) '_' num2str(SHbounds(2)) '.mat'], ...
     'data_B', 'V_Model_B', 'Model_B');

% Extract gravity anomaly in mGal
delta_g_mGal_B = data_B.vec.R * 1e5;


%% plot results
figure
aa = 18;
imagesc(longitudes, latitudes, delta_g_mGal);
c = colorbar;
ylabel(c, 'Gravity Anomaly (mGal)', 'Interpreter', 'latex', 'Fontsize', aa)
set(gca, 'YDir', 'normal', 'Fontsize', 12)
xlabel('Longitude ($^\circ$)', 'Interpreter', 'latex', 'Fontsize', aa)
ylabel('Latitude ($^\circ$)', 'Interpreter', 'latex', 'Fontsize', aa)
set(gca, 'ylim', [-90 90]);
set(gca, 'ytick', -90:30:90);
set(gca, 'xlim', [-180 180]);
set(gca, 'xtick', -180:30:180);


%% Compute degree variance for first model
% use function degree variance
% Extract the relevant columns
n_all = V_Model(:,1);  % Degree n
C_all = V_Model(:,3);  % Cnm
S_all = V_Model(:,4);  % Snm

nmax = max(n_all);
deg_var = zeros(nmax+1, 1); % Preallocate; deg 0 to nmax
degrees = 0:nmax;

% Loop over degrees and sum Cnm^2 + Snm^2 for each n
for n = degrees
    idx = (n_all == n);           % Find all rows with degree n
    deg_var(n+1) = sum(C_all(idx).^2 + S_all(idx).^2);  % Store in index n+1
end

% for degree power
deg_power = deg_var ./ (2*degrees + 1);
deg_power(1) = deg_var(1); % Avoid divide by zero for n = 0

%% Compute degree variance for Model B
% Extract the relevant columns
n_all_B = V_Model_B(:,1);  % Degree n
C_all_B = V_Model_B(:,3);  % Cnm
S_all_B = V_Model_B(:,4);  % Snm

nmax_B = max(n_all_B);
deg_var_B = zeros(nmax_B+1, 1); % Preallocate; deg 0 to nmax
degrees_B = 0:nmax_B;

% Loop over degrees and sum Cnm^2 + Snm^2 for each n
for n = degrees_B
    idx_B = (n_all_B == n);  % Find all rows with degree n
    deg_var_B(n+1) = sum(C_all_B(idx_B).^2 + S_all_B(idx_B).^2);  % Store in index n+1
end

% Compute degree power spectrum
deg_power_B = deg_var_B ./ (2*degrees_B + 1);
deg_power_B(1) = deg_var_B(1);  % Avoid divide by zero for n = 0


%% Plot the degree variance
figure;
semilogy(degrees, deg_power, 'b', 'LineWidth', 1.5); hold on;
semilogy(degrees_B, deg_power_B, 'r', 'LineWidth', 1.5);
xlabel('Spherical Harmonic Degree $n$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Degree Power', 'Interpreter', 'latex', 'FontSize', 14);
legend('Original Model', 'Model B', 'Location', 'northeast');
title('Degree Power Spectrum Comparison', 'FontSize', 16);
grid on;

%% Plot the degree variance
figure;
loglog(degrees, deg_var, 'b', 'LineWidth', 2); hold on;
loglog(degrees, deg_var_B, 'r', 'LineWidth', 2);
xlabel('Degree $n$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Degree Variance', 'Interpreter', 'latex', 'FontSize', 14);
legend();
title('Degree Variance Spectrum from Model SH Coefficients', 'FontSize', 14);
grid on;
