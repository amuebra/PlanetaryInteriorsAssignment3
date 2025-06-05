% Mercury Gravity Anomaly Map from MESSENGER SHADR file

% -------------------------------------------------------------------------
% PARAMETERS
% -------------------------------------------------------------------------
filename = 'ggmes_100v08_sha.tab';  % Path to SHA file
lmax = 100;                         % Maximum degree/order
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
latitudes = -90 : 1/resolution : 90;  % now from South to North
latitudes = latitudes(1:end-1);       % remove duplicate 90
longitudes = -180 : 1/resolution : 180;
longitudes = longitudes(1:end-1);     % remove duplicate 180

theta = deg2rad(90 - latitudes);      % colatitude in radians
phi = deg2rad(longitudes);            % longitude in radians
[phi_grid, theta_grid] = meshgrid(phi, theta); % (nlat x nlon)

% -------------------------------------------------------------------------
% GRAVITY ANOMALY COMPUTATION (\u2202V/\u2202r)
% -------------------------------------------------------------------------
delta_g = zeros(size(theta_grid));  % gravity anomaly in km/s^2

for l = 0:lmax
    Plm_all = zeros(l+1, length(theta));  % Precompute Plm for each theta
    for i = 1:length(theta)
        Plm_temp = legendre(l, cos(theta(i)), 'norm');
        Plm_all(:, i) = Plm_temp(:);
    end
    for m = 0:l
        P_lm = squeeze(Plm_all(m+1, :));
        P_grid = repmat(P_lm', 1, length(phi));

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
lat = latitudes;                 % Use directly, now -90 to 90
lon = rad2deg(phi);             % [1 x nlon]

% --- Shift longitudes so 180°E is at center (rotate planet)
lon_360 = mod(lon, 360);
[~, idx_180] = min(abs(lon_360 - 180));
lon_rotated = [lon_360(idx_180:end), lon_360(1:idx_180-1)];
delta_g_rotated = [delta_g(:, idx_180:end), delta_g(:, 1:idx_180-1)];
lon_wrapped = mod(lon_rotated + 180, 360) - 180;
[lon_grid, lat_grid] = meshgrid(lon_wrapped, lat);

fig = figure;
axesm('mollweid', ...
      'Frame', 'on', ...
      'Grid', 'on', ...
      'Origin', [0 0 0], ...
      'MapLatLimit', [-90 90], ...
      'MapLonLimit', [-180 180]);

pcolorm(lat_grid, lon_grid, delta_g_rotated);

c = colorbar;
c.Label.Interpreter = 'latex';
colormap(jet);
save('elevations.mat', 'elevations');
