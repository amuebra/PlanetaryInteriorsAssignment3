clear;
close all;
clc;

HOME = pwd;
addpath([HOME '/Data']);
addpath([HOME '/Results']);

% Iterative Bouguer Inversion to Determine Crustal Deviations (Δr)
% Requires gravity anomaly in gravity.m as `gravity_anomaly`

% Load gravity anomaly data
%run('Gravity.m');  % should define gravity_anomaly (in mGal)
load([HOME '/Results/gravity_anomaly_mGal.mat'], 'delta_g_mGal', 'latT', 'lonT')
longitudes = lonT;
latitudes = latT;

% === Parameters ===
D = 40e3;             % Average crustal thickness (m)
rho_c = 2800;         % Crustal density (kg/m^3)
rho_m = 3300;         % Mantle density (kg/m^3)
G = 6.67430e-11;      % Gravitational constant (m^3/kg/s^2)
scaling = 0.2;        % Scaling parameter (controls step size)
max_iter = 100;       % Max number of iterations
tolerance = 1e-3;     % Convergence criterion (in mGal)

% === Preprocess ===
% Convert mGal to m/s^2
obs_gravity = delta_g_mGal * 1e-5;
delta_rho = rho_m - rho_c;

% === Initial Model ===
sz = size(obs_gravity);
crust_thickness = D * ones(sz);  % Initial constant thickness

% === Iterative Inversion Loop ===
for iter = 1:max_iter
    % Compute gravity from current model
    % g = G * Δρ * Δr
    model_gravity = 2*pi*G * delta_rho * (crust_thickness - D);  % in m/s^2

    % Compute residual
    residual = obs_gravity - model_gravity;

    % Check convergence (use RMS)
    rms_residual = sqrt(mean(residual(:).^2)) / 1e-5;  % back to mGal
    fprintf('Iteration %d: RMS residual = %.4f mGal\n', iter, rms_residual);
    if rms_residual < tolerance
        disp('Converged!');
        break;
    end

    % Update model: Δr_new = Δr_old + scaling * residual / (G * Δρ)
    delta_r_update = scaling * residual / (2*pi*G * delta_rho);
    crust_thickness = crust_thickness + delta_r_update;

    % Optional: constrain model to physical range
    crust_thickness = max(crust_thickness, 0);  % no negative crust
end

% === Plot Final Model ===
figure;
imagesc(crust_thickness); colorbar;
title('Estimated Crustal Thickness (m)');
xlabel('Longitude'); ylabel('Latitude');
set(gca, 'YDir', 'normal');
% 
% % Optional: Save result
% save([HOME '/Results/crust_thickness_iterative.mat'], 'crust_thickness');

%% --- Plot
% fig = figure;
% axesm('mollweid', ...
%       'Frame', 'on', ...
%       'Grid', 'on', ...
%       'Origin', [0 0 0], ...   % CENTER AT 0°E!
%       'MapLatLimit', [-90 90], ...
%       'MapLonLimit', [-180 180]);
% 
% pcolorm(longitudes, latitudes, crust_thickness);
% 
% c = colorbar;
% c.Label.Interpreter = 'latex';
% %title('Mercury Gravity Anomaly (mGal) - Mollweide Projection');
% colormap(jet);