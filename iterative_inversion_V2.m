clear;
close all;
clc;

HOME = pwd;
addpath([HOME '/Data']);
addpath([HOME '/Results']);
addpath([HOME '/Tools'])

% Load the spherical harmonics coefficients from the file
filename = [HOME '/Data/ggmes_50v06_sha.tab'];
maxDegree = 20;
% Spherical harmonic synthesis settings
SHbounds = [1 maxDegree];
height = 0;
R_ref = 2439.4e3;       % Reference radius in meters
GM = 22031.815e9;       % Mercury GM (m^3/s^2)
G = 6.67430e-11;        % Gravitational constant (m^3/kg/s^2)
scaling = 1e6;
max_iter = 20;
tolerance = 1e-4;
coeffs = readmatrix(filename, 'FileType', 'text', 'Delimiter', ',');
coeffs = [[0,0,0,0,0,0]; coeffs];

% Load Topography
load([HOME '/Results/elevations.mat'], 'elevations')

% Use coefficients directly as V
V = coeffs;  % Already in [n m Cnm Snm] format
V(4,3) = 0; %set C20 to zero 

% Optional: limit to maxDegree
V = V(V(:,1) <= maxDegree, :);

% Define grid resolution
resolution = 1;
latLimT = [-90+(1/resolution/2) 90-(1/resolution/2) 1/resolution]; 
lonLimT = [1/resolution/2 360-(1/resolution/2) 1/resolution];     % 0 to 360 degree
%lonLimT = [-180+(1/resolution/2) 180-(1/resolution/2) 1/resolution]; %-180 to 180 degree

% Generate lat/lon grid
lonT = lonLimT(1):lonLimT(3):lonLimT(2);
%latT = fliplr(latLimT(1):latLimT(3):latLimT(2));
latT = latLimT(1):latLimT(3):latLimT(2);

%% computate model
% Construct new model
Model = struct();
Model.number_of_layers = 2;
Model.name = 'Mercury';

% Additional variables
Model.GM = GM;
Model.Re = R_ref;
Model.geoid = 'none';
Model.nmax = maxDegree;   
Model.correct_depth = 0;
D = 16000;

% Top layer (Crust)
Model.l1.bound = elevations;    % meters with respect to reference sphere
Model.l1.dens  = 2800;          % Density in kg/m3

% Second layer (Mantle)
Model.l2.bound = -D*ones(size(elevations)); % meters with respect to reference sphere
Model.l2.dens  = 3200;	        % Density in kg/m3

% Bottom bound
Model.l3.bound = -100000;       % meters with respect to reference sphere

% Global Spherical Harmonic Analysis 
[V_Model] = segment_2layer_model(Model.l1.bound,Model.l2.bound,Model.l3.bound,Model.l1.dens,Model.l2.dens,25000,Model);
V_Model(1,3) = 0;
V_Model(3,3) = 0;

% Run synthesis
model_result = model_SH_synthesis(lonLimT, latLimT, height, SHbounds, V_Model, Model);

% Extract gravity anomaly (in mGal)
deltag_model = model_result.vec.R;
deltag_model_mGal = model_result.vec.R * 1e5;

%% computate observations
% Define Model structure
Observation.GM = GM;
Observation.Re = R_ref;

% Run synthesis
observation_result = model_SH_synthesis(lonLimT, latLimT, height, SHbounds, V, Observation);

% Extract gravity anomaly (in mGal)
deltag_observation = observation_result.vec.R;
deltag_observation_mGal = observation_result.vec.R * 1e5;

%Constants
delta_rho = Model.l2.dens - Model.l1.dens;
residual_history = [];

%% === Iterative Inversion Loop ===
for iter = 1:max_iter
    % Compute gravity from current model
    [V_Model] = segment_2layer_model(Model.l1.bound,Model.l2.bound,Model.l3.bound,Model.l1.dens,Model.l2.dens,25000,Model);
    V_Model(1,3) = 0;
    V_Model(3,3) = 0;
    model_result = model_SH_synthesis(lonLimT, latLimT, height, SHbounds, V_Model, Model);
    deltag_model = model_result.vec.R;

    % compute residual
    residual = deltag_observation - flipud(deltag_model);
    disp(residual);
    
    % Check convergence (use RMS)
    rms_residual = max(residual(:))  % back to mGal
    residual_history(end+1) = rms_residual;
    fprintf('Iteration %d: RMS residual = %.4f m\s^2', iter, rms_residual);
    if abs(max(residual(:))) < tolerance
         disp('Converged!');
         break;
    end

    % update model
    %delta_r_update = scaling * residual/(2*pi*G*delta_rho);
    delta_r_update = scaling.*residual;
    Model.l2.bound = Model.l2.bound + delta_r_update;
        
   
    % figure;
    % aa = 18;
    % imagesc(lonT, latT, deltag_model_mGal);
    % c = colorbar;
    % ylabel(c, 'Gravity Anomaly (mGal)', 'Interpreter', 'latex', 'Fontsize', aa)
    % set(gca, 'YDir', 'normal', 'Fontsize', 12)
    % xlabel('Longitude ($^\circ$)', 'Interpreter', 'latex', 'Fontsize', aa)
    % ylabel('Latitude ($^\circ$)', 'Interpreter', 'latex', 'Fontsize', aa)
    % set(gca, 'ylim', [-90 90]);
    % set(gca, 'ytick', -90:30:90);
    % set(gca, 'xlim', [-180 180]);
    % set(gca, 'xtick', -180:30:180);
end

%% plot results
figure;
plot(1:length(residual_history), residual_history, 'o-', 'LineWidth', 2);
xlabel('Iteration Number', 'FontSize', 12, 'Interpreter', 'latex');
ylabel('RMS Residual (mGal)', 'FontSize', 12, 'Interpreter', 'latex');
title('Residual vs. Iteration', 'FontSize', 14, 'Interpreter', 'latex');
grid on;
figure;
aa = 18;
imagesc(lonT, latT, flipud(deltag_model)./1e-5);
c = colorbar;
ylabel(c, 'Gravity Anomaly Model(mGal)', 'Interpreter', 'latex', 'Fontsize', aa)
set(gca, 'YDir', 'normal', 'Fontsize', 12)
xlabel('Longitude ($^\circ$)', 'Interpreter', 'latex', 'Fontsize', aa)
ylabel('Latitude ($^\circ$)', 'Interpreter', 'latex', 'Fontsize', aa)
set(gca, 'ylim', [-90 90]);
set(gca, 'ytick', -90:30:90);
%set(gca, 'xlim', [-180 180]);
%set(gca, 'xtick', -180:30:180);

figure;
aa = 18;
imagesc(lonT, latT, deltag_observation./1e-5);
c = colorbar;
ylabel(c, 'Gravity Anomaly Observation (mGal)', 'Interpreter', 'latex', 'Fontsize', aa)
set(gca, 'YDir', 'normal', 'Fontsize', 12)
xlabel('Longitude ($^\circ$)', 'Interpreter', 'latex', 'Fontsize', aa)
ylabel('Latitude ($^\circ$)', 'Interpreter', 'latex', 'Fontsize', aa)
set(gca, 'ylim', [-90 90]);
set(gca, 'ytick', -90:30:90);
%set(gca, 'xlim', [-180 180]);
%set(gca, 'xtick', -180:30:180);

figure;
aa = 18;
imagesc(lonT, latT, residual./1e-5);
c = colorbar;
ylabel(c, 'residual (mGal)', 'Interpreter', 'latex', 'Fontsize', aa)
set(gca, 'YDir', 'normal', 'Fontsize', 12)
xlabel('Longitude ($^\circ$)', 'Interpreter', 'latex', 'Fontsize', aa)
ylabel('Latitude ($^\circ$)', 'Interpreter', 'latex', 'Fontsize', aa)
set(gca, 'ylim', [-90 90]);
set(gca, 'ytick', -90:30:90);
%%
figure;
aa = 18;
thickness = elevations - Model.l2.bound;
fprintf('Minimum thickness = %.4f meters\n', min(thickness(:)));
imagesc(lonT, latT, thickness);
c = colorbar;
ylabel(c, 'L2 Bound (m)', 'Interpreter', 'latex', 'Fontsize', aa)
set(gca, 'YDir', 'normal', 'Fontsize', 12)
xlabel('Longitude ($^\circ$)', 'Interpreter', 'latex', 'Fontsize', aa)
ylabel('Latitude ($^\circ$)', 'Interpreter', 'latex', 'Fontsize', aa)
set(gca, 'ylim', [-90 90]);
set(gca, 'ytick', -90:30:90);
colormap(turbo);
%set(gca, 'xlim', [-180 180]);
%set(gca, 'xtick', -180:30:180);

%%
figure;
aa = 18;
imagesc(lonT, latT, elevations);
c = colorbar;
ylabel(c, 'Elevation', 'Interpreter', 'latex', 'Fontsize', aa)
set(gca, 'YDir', 'normal', 'Fontsize', 12)
xlabel('Longitude ($^\circ$)', 'Interpreter', 'latex', 'Fontsize', aa)
ylabel('Latitude ($^\circ$)', 'Interpreter', 'latex', 'Fontsize', aa)
set(gca, 'ylim', [-90 90]);
set(gca, 'ytick', -90:30:90);
colormap(turbo);
%% redo Model
% % update model
% Model.l2.bound = Model.l2.bound + delta_r_update;
% 
% % Global Spherical Harmonic Analysis 
% [V_Model] = segment_2layer_model(Model.l1.bound,Model.l2.bound,Model.l3.bound,Model.l1.dens,Model.l2.dens,25000,Model);
% 
% % Run synthesis
% model_result = model_SH_synthesis(lonLimT, latLimT, height, SHbounds, V_Model, Model);
% 
% % Extract gravity anomaly (in mGal)
% deltag_model_mGal = model_result.vec.R * 1e5;

%% plot results
% figure;
% aa = 18;
% imagesc(lonT, latT, deltag_model_mGal);
% c = colorbar;
% ylabel(c, 'Gravity Anomaly (mGal)', 'Interpreter', 'latex', 'Fontsize', aa)
% set(gca, 'YDir', 'normal', 'Fontsize', 12)
% xlabel('Longitude ($^\circ$)', 'Interpreter', 'latex', 'Fontsize', aa)
% ylabel('Latitude ($^\circ$)', 'Interpreter', 'latex', 'Fontsize', aa)
% set(gca, 'ylim', [-90 90]);
% set(gca, 'ytick', -90:30:90);
% set(gca, 'xlim', [-180 180]);
% set(gca, 'xtick', -180:30:180);

%% save data
%save([HOME '/Results/deltag_model_mGal.mat'], 'deltag_model_mGal', 'latT', 'lonT');
%save([HOME '/Results/deltag_observagtion_mGal.mat'], 'deltag_observation_mGal', 'latT', 'lonT');