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
G = 6.67430e-11;        % Gravitational constant (m^3/kg/s^2)
scaling = -0.01;
max_iter = 1;
coeffs = readmatrix(filename, 'FileType', 'text', 'Delimiter', ',');

% Load Topography
load([HOME '/Results/elevations.mat'], 'elevations')

% Use coefficients directly as V
V = coeffs;  % Already in [n m Cnm Snm] format
V(3,3) = 0; %set C20 to zero 

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

%% computate model
% Construct new model
Model = struct();
Model.number_of_layers = 2;
Model.name = 'Mercury';

% Additional variables
Model.GM = GM;
Model.Re = R_ref;
Model.geoid = 'none';
Model.nmax = 20;     
Model.correct_depth = 0;

% Top layer (Crust)
Model.l1.bound = elevations;    % meters with respect to reference sphere
Model.l1.dens  = 2850;          % Density in kg/m3

% Second layer (Mantle)
Model.l2.bound = -29360*ones(size(elevations)); % meters with respect to reference sphere
Model.l2.dens  = 3100;	        % Density in kg/m3

% Bottom bound
Model.l3.bound = -100000;       % meters with respect to reference sphere

% Global Spherical Harmonic Analysis 
[V_Model] = segment_2layer_model(Model.l1.bound,Model.l2.bound,Model.l3.bound,Model.l1.dens,Model.l2.dens,25000,Model);
V_Model(1,3) = 0;
V_Model(3,3) = 0;
% Spherical harmonic synthesis settings
SHbounds = [1 Model.nmax];
height = 0;

% Run synthesis
model_result = model_SH_synthesis(lonLimT, latLimT, height, SHbounds, V_Model, Model);

% Extract gravity anomaly (in mGal)
deltag_model_mGal = model_result.vec.R * 1e5;

%% computate observations
% Define Model structure
Observation.GM = GM;
Observation.Re = R_ref;

% Spherical harmonic synthesis settings
SHbounds = [1 20];
height = 0;

% Run synthesis
observation_result = model_SH_synthesis(lonLimT, latLimT, height, SHbounds, V, Observation);

% Extract gravity anomaly (in mGal)
deltag_observagtion_mGal = observation_result.vec.R * 1e5;

%% === Iterative Inversion Loop ===
for iter = 1:max_iter
    % compute residual
    residual_mGal = deltag_observagtion_mGal - deltag_model_mGal;
    delta_rho = Model.l2.dens - Model.l1.dens;
    delta_r_update = scaling * residual_mGal * 1e-5 / (2*pi*G * delta_rho);
    % delta_r_update = residual_mGal * 1e-5 / (2*pi*G * delta_rho);
    
    % Check convergence (use RMS)
    rms_residual = sqrt(mean(residual_mGal(:).^2)) / 1e-5;  % back to mGal
    fprintf('Iteration %d: RMS residual = %.4f mGal\n', iter, rms_residual);
    % if rms_residual < tolerance
    %     disp('Converged!');
    %     break;
    % end

    % update model
    Model.l2.bound = Model.l2.bound - delta_r_update;
    
    % Global Spherical Harmonic Analysis 
    [V_Model] = segment_2layer_model(Model.l1.bound,Model.l2.bound,Model.l3.bound,Model.l1.dens,Model.l2.dens,25000,Model);
    V_Model(1,3) = 0;
    V_Model(3,3) = 0;

    % Run synthesis
    model_result = model_SH_synthesis(lonLimT, latLimT, height, SHbounds, V_Model, Model);
    
    % Extract gravity anomaly (in mGal)
    deltag_model_mGal = model_result.vec.R * 1e5;
    
    figure;
    aa = 18;
    imagesc(lonT, latT, deltag_model_mGal);
    c = colorbar;
    ylabel(c, 'Gravity Anomaly (mGal)', 'Interpreter', 'latex', 'Fontsize', aa)
    set(gca, 'YDir', 'normal', 'Fontsize', 12)
    xlabel('Longitude ($^\circ$)', 'Interpreter', 'latex', 'Fontsize', aa)
    ylabel('Latitude ($^\circ$)', 'Interpreter', 'latex', 'Fontsize', aa)
    set(gca, 'ylim', [-90 90]);
    set(gca, 'ytick', -90:30:90);
    set(gca, 'xlim', [-180 180]);
    set(gca, 'xtick', -180:30:180);
end

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
% 
% figure;
% aa = 18;
% imagesc(lonT, latT, deltag_observagtion_mGal);
% c = colorbar;
% ylabel(c, 'Gravity Anomaly (mGal)', 'Interpreter', 'latex', 'Fontsize', aa)
% set(gca, 'YDir', 'normal', 'Fontsize', 12)
% xlabel('Longitude ($^\circ$)', 'Interpreter', 'latex', 'Fontsize', aa)
% ylabel('Latitude ($^\circ$)', 'Interpreter', 'latex', 'Fontsize', aa)
% set(gca, 'ylim', [-90 90]);
% set(gca, 'ytick', -90:30:90);
% set(gca, 'xlim', [-180 180]);
% set(gca, 'xtick', -180:30:180);
% 
% figure;
% aa = 18;
% imagesc(lonT, latT, residual_mGal);
% c = colorbar;
% ylabel(c, 'residual (mGal)', 'Interpreter', 'latex', 'Fontsize', aa)
% set(gca, 'YDir', 'normal', 'Fontsize', 12)
% xlabel('Longitude ($^\circ$)', 'Interpreter', 'latex', 'Fontsize', aa)
% ylabel('Latitude ($^\circ$)', 'Interpreter', 'latex', 'Fontsize', aa)
% set(gca, 'ylim', [-90 90]);
% set(gca, 'ytick', -90:30:90);
% set(gca, 'xlim', [-180 180]);
% set(gca, 'xtick', -180:30:180);

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
save([HOME '/Results/deltag_model_mGal.mat'], 'deltag_model_mGal', 'latT', 'lonT');
save([HOME '/Results/deltag_observagtion_mGal.mat'], 'deltag_observagtion_mGal', 'latT', 'lonT');