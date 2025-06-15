%% Initialization
clear; close all; clc;

HOME = pwd;
addpath([HOME '/Data']);
addpath([HOME '/Results']);
addpath([HOME '/Tools']);
addpath([HOME '/ScientificColourMaps8/vik']);
load('vik.mat'); % Load colormap

%% Load Data
load([HOME '/Results/Airy_thickness.mat'], 'crust_thickness');
load([HOME '/Results/coeffs_obs.mat'], 'V');
load([HOME '/Results/elevations.mat'], 'elevations');

% Model constants
resolution = 1;
height = 0;
Model.GM = 22031.815E9;      % Gravitational constant
Model.Re = 2439.4E3;         % Reference radius
Model.geoid = 'none';
Model.nmax = 50;
Model.correct_depth = 0;
SHbounds = [1 50];
nu = 0.25;  % Poisson’s ratio

% Grid setup
latLim = [-90 + (1/resolution/2), 90 - (1/resolution/2), 1/resolution];
lonLim = [0 + (1/resolution/2), 360 - (1/resolution/2), 1/resolution];

lonT = lonLim(1):lonLim(3):lonLim(2);
latT = latLim(1):latLim(3):latLim(2);
LonT = repmat(lonT, length(latT), 1);
LatT = repmat(latT', 1, length(lonT));

aa = 16; % Font size

%% GSHA (Airy model coefficients)
cs = GSHA(crust_thickness, 179);  % Spherical harmonic analysis
sc = cs2sc(cs);                   % Convert to SH coefficients
n = 1:size(sc, 1);                % Degrees

%% Flexural Response for Different Te Values
Te_values = [10e3, 20e3, 30e3, 40e3, 50e3];  % Elastic thicknesses to test [m]
flex_maps = zeros(length(latT), length(lonT), length(Te_values));  % Store results

for k = 1:length(Te_values)
    Te = Te_values(k);
    
    % Flexural rigidity
    D = 100e9 * Te^3 / (12 * (1 - nu^2));
    
    % Flexural filter
    PHI = (1 + (D / (400 * 3.7)) .* (2 * (n + 1) / (2 * Model.Re)).^4) .^ (-1);
    R = 2439.4e3;
    %PHI = (1+D/(400*3.7).*(1/R.^4.*(n.*(n+1)-2).^2./(1-(1-nu.^2)./(n.*(n+1)))+12.*(1-nu)./(Te.^2.*R.^2)*(1-2./(n.*(n+1))./(1-(1-nu)./(n.*(n+1)))))).^(-1);

    
    % Apply filter to SH coefficients
    sc_flex = zeros(size(sc));
    for m = 1:size(sc, 2)
        sc_flex(:, m) = sc(:, m) .* PHI';
    end
    
    % Synthesize spatial map
    flex_maps(:,:,k) = GSHS(sc_flex, lonT, 90 - latT, 179);
end

%% Plot Comparison of Flexural Thickness for Each Te
figure('Position', [100 100 1600 300]);
for k = 1:length(Te_values)
    subplot(1, length(Te_values), k)
    imagesc(lonT, latT, flex_maps(:,:,k) ./ 1e3);  % Convert to km
    set(gca, 'YDir', 'normal', 'FontSize', aa);
    xlabel('Longitude ($^\circ$)', 'Interpreter', 'latex', 'FontSize', aa);
    if k == 1
        ylabel('Latitude ($^\circ$)', 'Interpreter', 'latex', 'FontSize', aa);
    end
    title(['T$_e$ = ' num2str(Te_values(k)/1e3) ' km'], 'Interpreter', 'latex', 'FontSize', aa);
    colormap(vik);
    axis equal tight;
    set(gca, 'ylim', [-90 90], 'ytick', -90:30:90);
    set(gca, 'xlim', [0 360], 'xtick', 0:60:360);
end

% Shared colorbar
hp4 = get(subplot(1,length(Te_values),length(Te_values)),'Position');
h = colorbar('Position', [hp4(1)+hp4(3)+0.01, hp4(2), 0.015, hp4(2)+hp4(3)*1.5]);
ylabel(h, 'Flexural Crust Thickness (km)', 'Interpreter', 'latex', 'FontSize', aa);


%% Compute Degree Variance for Multiple Te Values

% Observation data
V = sortrows(V, [2, 1]);
[n_obs, DV_obs] = degreeVariance(V);

DV_models = cell(length(Te_values), 1); % Store DV for each Te
n_models = cell(length(Te_values), 1);  % Store degrees

for k = 1:length(Te_values)
    % Setup model structure
    Model_B = Model;  % Start from base model
    Model_B.name = ['Flexural_Te_' num2str(Te_values(k)/1e3) 'km'];
    
    % Model layers
    Model_B.l1.bound = elevations;
    Model_B.l1.dens  = 2800;
    Model_B.l2.dens  = 3200;
    Model_B.l3.bound = -100000;
    Model_B.l2.bound = -flex_maps(:,:,k);  % Set crust-mantle boundary

    % SH coefficients
    [V_Model_B] = segment_2layer_model(Model_B.l1.bound, Model_B.l2.bound, Model_B.l3.bound, ...
                                       Model_B.l1.dens, Model_B.l2.dens, 62000, Model_B);
    % Remove degree-0 and degree-2 zonal terms (optional)
    V_Model_B(1,3) = 0;  % n=0
    V_Model_B(3,3) = 0;  % n=2, m=0

    % Compute degree variance
    [n_B, DV_B] = degreeVariance(V_Model_B);
    DV_models{k} = DV_B;
    n_models{k} = n_B;
end

%% Comparison between curves
function [rmse_values] = compute_dv_misfit(n_obs, DV_obs, n_models, DV_models)
% Compute RMSE between observed and model degree variance (in log space)

rmse_values = zeros(length(DV_models), 1);

% Remove invalid obs (NaNs or zero/negative)
valid_obs_idx = ~isnan(DV_obs) & DV_obs > 0;
n_obs_valid = n_obs(valid_obs_idx);
DV_obs_valid = DV_obs(valid_obs_idx);

for k = 1:length(DV_models)
    n_mod = n_models{k};
    DV_mod = DV_models{k};

    % Remove invalid model entries
    valid_mod_idx = ~isnan(DV_mod) & DV_mod > 0;
    n_mod_valid = n_mod(valid_mod_idx);
    DV_mod_valid = DV_mod(valid_mod_idx);

    % Interpolate model to obs degrees
    DV_mod_interp = interp1(n_mod_valid, DV_mod_valid, n_obs_valid, 'linear', 'extrap');

    % Remove any new NaNs or invalid interpolated values
    valid_interp = ~isnan(DV_mod_interp) & DV_mod_interp > 0;

    if any(valid_interp)
        diff = log10(DV_mod_interp(valid_interp)) - log10(DV_obs_valid(valid_interp));
        rmse_values(k) = sqrt(mean(diff.^2));
    else
        rmse_values(k) = NaN;
    end
end

end


%% Compute RMSE between models and observation
degree_mask = (n_obs >= 6 & n_obs <= 50);
% Filter the observations
n_obs_masked = n_obs(degree_mask);
DV_obs_masked = DV_obs(degree_mask);

% Then call the function with filtered data
rmse_values = compute_dv_misfit(n_obs_masked, DV_obs_masked, n_models, DV_models);

% Print results
fprintf('\nRMSE (log10 scale) between models and observation:\n');
for k = 1:length(Te_values)
    fprintf('T_e = %2d km: RMSE = %.4f\n', Te_values(k)/1e3, rmse_values(k));
end

% Find best-fitting Te
[~, idx_best] = min(rmse_values);
best_Te = Te_values(idx_best);
fprintf('\n>> Best-fit T_e = %d km (minimum RMSE)\n', best_Te/1e3);

% Define degree range of interest
n_min = 6;
n_max = 50;

% Filter observation data
mask_obs = (n_obs >= n_min) & (n_obs <= n_max) & ~isnan(DV_obs) & (DV_obs > 0);
n_obs_filt = n_obs(mask_obs);
DV_obs_filt = DV_obs(mask_obs);

%% Start plotting
figure('units', 'points', 'Position', [0, 0, 600, 400]);

% Color setup
colors = lines(length(Te_values));
hold on;

% Plot each model
for k = 1:length(Te_values)
    n_mod = n_models{k};
    DV_mod = DV_models{k};

    % Filter model data
    mask_mod = (n_mod >= n_min) & (n_mod <= n_max) & ~isnan(DV_mod) & (DV_mod > 0);
    n_mod_filt = n_mod(mask_mod);
    DV_mod_filt = DV_mod(mask_mod);

    % Optional: Interpolate model to obs degrees (to match markers)
    % DV_mod_interp = interp1(n_mod_filt, DV_mod_filt, n_obs_filt, 'linear', 'extrap');

    scatter(n_mod_filt, DV_mod_filt, 36, 'MarkerEdgeColor', colors(k,:), ...
        'DisplayName', ['T$_e$ = ' num2str(Te_values(k)/1e3) ' km'], ...
        'LineWidth', 1.2);
end

% Plot filtered observations
scatter(n_obs_filt, DV_obs_filt, 50, 'k^', 'LineWidth', 1.5, 'DisplayName', 'Observation');

% Axis formatting
xlabel('Degree $n$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Degree Variance', 'Interpreter', 'latex', 'FontSize', 14);
set(gca, 'YScale', 'log', 'XScale', 'log');
legend('Location', 'southwest', 'Interpreter', 'latex');
grid on;
% title('Degree Variance Comparison for Different $T_e$', 'FontSize', 14, 'Interpreter', 'latex');

% Save figure
saveas(gcf, 'Figures/FlexuralTeChange.svg');




