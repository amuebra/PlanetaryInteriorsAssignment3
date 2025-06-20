% main file for the complete GSH circle for a particular model
clear;
close all;
clc;


HOME = pwd;
addpath([HOME '/Data']);
addpath([HOME '/Results'])
addpath([HOME '/Tools']);

% Load Data

load([HOME '/Results/elevations.mat'], 'elevations');
load([HOME '/Results/Airy_thickness.mat'], 'crust_thickness','longitudes', 'latitudes');
airy_thickness = crust_thickness;
load([HOME '/Results/Flexural_thickness.mat'], 'mapf');
flexural_crust_thickness = mapf;
load([HOME '/Results/coeffs_obs.mat'], 'V');
load([HOME '/Results/coeffs_Bouguer.mat'], 'V_Model_C');

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
  Model.nmax = 50;     
  Model.correct_depth = 0;
  
  % Top layer
  %Topography (maxtrix latitude and logitude in m) 
  %Model.l1.bound = gmt2matrix(load([HOME '/Data/crust1.bd1.gmt'])).*1e3;  % meters with respect to reference sphere
  Model.l1.bound = elevations;
  Model.l1.dens  = 2800;
  
  %% Second layer
  Model.l2.bound = -airy_thickness;
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
lonLim =    [0.5 359.5 1];% [deg] min longitude, max longitude, resolution longitude (preferable similar to latitude)
height =    0; % height of computation above spheroid
SHbounds =  [0 50]; % Truncation settings: lower limit, upper limit SH-coefficients used

%%%%%%%%%%%%%% Part that can be modified %%%%%%%%%%%%%%%%%%%%%%%

%% Global Spherical Harmonic Analysis 

tic;
[V_Model] = segment_2layer_model(Model.l1.bound,Model.l2.bound,Model.l3.bound,Model.l1.dens,Model.l2.dens,62000,Model); % order, number, S-coefficentes, c-coefficents
V_Model(1,3) = 0;
V_Model(3,3)= 0;
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
V_Model_B(1,3) = 0;
V_Model_B(3,3) = 0;

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
%set(gca, 'xlim', [-180 180]);
%set(gca, 'xtick', -180:30:180);

%% Compute degree variance
[n_A, DV_ModelA] = degreeVariance(V_Model);
[n_B, DV_ModelB] = degreeVariance(V_Model_B);
[n_C, DV_ModelC] = degreeVariance(V_Model_C);
V = sortrows(V, [2, 1]);
[n, DV] = degreeVariance(V);

%% Plot the degree variance
figure('units', 'points', 'Position', [0, 0, 455.2441, 0.5*455.2441]);
scatter(n_A.', DV_ModelA', 'b', 'LineWidth', 2); hold on;
scatter(n_B.', DV_ModelB', 'r', 'LineWidth', 2);
scatter(n_C.', DV_ModelC', 'g', 'LineWidth', 2);
scatter(n.', DV', 'k', '^', 'LineWidth', 2);
xlabel('Degree $n$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Degree Variance', 'Interpreter', 'latex', 'FontSize', 14);
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
legend('Airy Model', 'Flexural Model', 'Bouguer Model', 'Observation');
% title('Degree Variance Spectrum from Model SH Coefficients', 'FontSize', 14);
grid on;
saveas(gcf, 'Figures/Degree_Variance.svg');

%% Compute degree power spectrum
deg_power_A = DV_ModelA ./ (2*n_A + 1);
deg_power_B = DV_ModelB ./ (2*n_B + 1);
deg_power_C = DV_ModelC ./ (2*n_C + 1);
deg_power = DV ./ (2*n + 1);
% deg_power_A_mGal = deg_power_A*1e5;
% deg_power_B_mGal = deg_power_B*1e5;

figure('units', 'points', 'Position', [0, 0, 455.2441, 0.5*455.2441]);
scatter(n_A.', deg_power_A', 'b', 'LineWidth', 2); hold on;
scatter(n_B.', deg_power_B', 'r', 'LineWidth', 2);
scatter(n_C.', deg_power_C', 'g', 'LineWidth', 2);
scatter(n.', deg_power', 'k', '^', 'LineWidth', 2);
xlabel('Degree $n$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Power Spectrum', 'Interpreter', 'latex', 'FontSize', 14);
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
legend('Airy Model', 'Flexural Model', 'Bouguer Model', 'Observation');
%title('Degree Variance Spectrum from Model SH Coefficients', 'FontSize', 14);
grid on;
saveas(gcf, 'Figures/Power_Spectrum.svg');

% figure;
% semilogy(degrees, deg_power, 'b', 'LineWidth', 1.5); hold on;
% semilogy(degrees_B, deg_power_B, 'r', 'LineWidth', 1.5);
% xlabel('Spherical Harmonic Degree $n$', 'Interpreter', 'latex', 'FontSize', 14);
% ylabel('Degree Power', 'Interpreter', 'latex', 'FontSize', 14);
% legend('Original Model', 'Model B', 'Location', 'northeast');
% title('Degree Power Spectrum Comparison', 'FontSize', 16);
% grid on;

%% Compute Root Mean Square
RMS_A = sqrt(deg_power_A);
RMS_B = sqrt(deg_power_B);
RMS_C = sqrt(deg_power_C);
RMS = sqrt(deg_power);

figure('units', 'points', 'Position', [0, 0, 455.2441, 0.5*455.2441]);
scatter(n_A.', RMS_A', 'b', 'LineWidth', 2); hold on;
scatter(n_B.', RMS_B', 'r', 'LineWidth', 2);
scatter(n_C.', RMS_C', 'g', 'LineWidth', 2);
scatter(n.', RMS', 'k', '^', 'LineWidth', 2);
xlabel('Degree $n$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Root Mean Square', 'Interpreter', 'latex', 'FontSize', 14);
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
legend('Airy Model', 'Flexural Model', 'Bouguer Model', 'Observation');
%title('Degree Variance Spectrum from Model SH Coefficients', 'FontSize', 14);
grid on;
saveas(gcf, 'Figures/Root_Mean_Square.svg');