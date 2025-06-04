% Mercury Gravity Anomaly Map from MESSENGER SHADR file

% -------------------------------------------------------------------------
% PARAMETERS
% -------------------------------------------------------------------------
filename = 'ggmes_20v04_sha.tab';  % Path to SHA file
lmax = 20;                         % Maximum degree/order
R_ref = 2439.4;                     % Reference radius (km)
GM = 22031.8150000000;              % Mercury GM (km^3/s^2)
resolution = 1;                     % degrees (1 = 1x1°, 4 = 0.25°)

% -------------------------------------------------------------------------
% READ COEFFICIENTS
% -------------------------------------------------------------------------
fid = fopen(filename, 'r');
if fid == -1
    error('Cannot open file: %s', filename);
end

% Skip header (first 2 lines)
fgetl(fid); fgetl(fid);

% Read coefficients: degree, order, C, S, C_unc, S_unc
data = textscan(fid, '%f%f%f%f%f%f', 'Delimiter', ',', 'CollectOutput', true);
fclose(fid);

coeffs = data{1};
degree = coeffs(:,1);
order  = coeffs(:,2);
Clm    = coeffs(:,3);
Slm    = coeffs(:,4);

% Initialize coefficient matrices
Clm_mat = zeros(lmax+1);
Slm_mat = zeros(lmax+1);
for i = 1:length(degree)
    l = degree(i);
    m = order(i);
    if l <= lmax && m <= l
        Clm_mat(l+1, m+1) = Clm(i);
        Slm_mat(l+1, m+1) = Slm(i);
    end
end

% -------------------------------------------------------------------------
% LAT/LON GRID (user-defined resolution)
% -------------------------------------------------------------------------
latitudes = 90 : -1/resolution : -90;  % now from North to South
latitudes = latitudes(1:end-1);  % remove duplicate 90
longitudes = -180 : 1/resolution : 180;
longitudes = longitudes(1:end-1); % remove duplicate 180

theta = deg2rad(90 - latitudes);      % colatitude in radians
phi = deg2rad(longitudes);            % longitude in radians
[phi_grid, theta_grid] = meshgrid(phi, theta); % (nlat x nlon)

% -------------------------------------------------------------------------
% GRAVITY ANOMALY COMPUTATION (∂V/∂r)
% -------------------------------------------------------------------------
delta_g = zeros(size(theta_grid));  % gravity anomaly in km/s^2

for l = 0:lmax
    Plm_all = legendre(l, cos(theta), 'norm');  % normalized, (m+1 x nlat)
    for m = 0:l
        P_lm = squeeze(Plm_all(m+1, :));         % 1 x nlat
        P_grid = repmat(P_lm', 1, length(phi));  % nlat x nlon

        delta_g = delta_g + (l + 1)* (R_ref / R_ref)^(l + 2) * P_grid .* ...
            (Clm_mat(l+1, m+1) * cos(m * phi_grid) + ...
             Slm_mat(l+1, m+1) * sin(m * phi_grid));
    end
end

% Scale and convert to mGal
delta_g_kms2 = GM / R_ref^2 * delta_g;
delta_g_mGal = delta_g_kms2 * 1e8;  % 1 km/s^2 = 1e8 mGal

% -------------------------------------------------------------------------
% PLOT GRAVITY MAP
% -------------------------------------------------------------------------
 % figure;
 % imagesc(longitudes, latitudes, delta_g_mGal);
 % axis xy;
 % colorbar;
 % title('Mercury Gravity Anomaly from MESSENGER (mGal)');
 % xlabel('Longitude (°)');
 % ylabel('Latitude (°)');
 % colormap jet;

% % Assume you have delta_g [nlat x nlon], theta [1 x nlat], phi [1 x nlon]
% 
% Convert colatitude (theta) to latitude
lat = 90 - rad2deg(theta);  % latitude in degrees
lon = rad2deg(phi);         % longitude in degrees

% Ensure longitude is 0–360, then wrap to -180 to 180 for plotting
lon_wrapped = mod(lon + 180, 360) - 180;

% Sort longitudes so that they are in increasing order (-180 to 180)
[lon_sorted, idx] = sort(lon_wrapped);
delta_g_sorted = delta_g(:, idx);

% Create meshgrid for mapping
[lon_grid, lat_grid] = meshgrid(lon_sorted, lat);

% Plot using Mapping Toolbox
figure;
axesm ('mollweid', 'Frame', 'on', 'Grid', 'on', ...
       'Origin', [0 180 0], ...      % Central meridian at 180°E
       'MapLatLimit', [-90 90], ...
       'MapLonLimit', [-180 180]);

% Use surfm or pcolorm to plot
pcolorm(lat_grid, lon_grid, delta_g_sorted);

% Add coastlines if needed
%coast = load('coastlines');
%plotm(coast.coastlat, coast.coastlon, 'k');

% Add colorbar and title
colorbar;
title('Gravity Anomaly (m/s^2 or mGal) - Mollweide Projection');
