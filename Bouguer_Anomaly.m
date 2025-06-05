% Mercury Gravity Anomaly Map from MESSENGER SHADR file

% -------------------------------------------------------------------------
%% PARAMETERS
% -------------------------------------------------------------------------
filename = 'ggmes_100v08_sha.tab';  % Path to SHA file
lmax = 100;                         % Maximum degree/order
R_ref = 2439.4;                     % Reference radius (km)
GM = 22031.8150000000;              % Mercury GM (km^3/s^2)
resolution = 1;                     % degrees (1 = 1x1°, 4 = 0.25x0.25°)

% -------------------------------------------------------------------------
%% READ COEFFICIENTS
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
%% READ TOPOGRAPHY DATA
% -------------------------------------------------------------------------
data = read(Tiff('Mercury_Messenger_USGS_DEM_Global_665m_v2.tif'));

resolution = 1;
latitudes = -90 : 1/resolution : 90;
latitudes = latitudes(1:end-1);
longitudes = -180 : 1/resolution : 180;
longitudes = longitudes(1:end-1);

multiplier = 0.5;
elevations = multiplier * double(data);

% -------------------------------------------------------------------------
%% LAT/LON GRID (user-defined resolution)
% -------------------------------------------------------------------------
latitudes = 90 : -1/resolution : -90;  % now from North to South
latitudes = latitudes(1:end-1);  % remove duplicate 90
longitudes = -180 : 1/resolution : 180;
longitudes = longitudes(1:end-1); % remove duplicate 180

theta = deg2rad(90 - latitudes);      % colatitude in radians
phi = deg2rad(longitudes);            % longitude in radians
[phi_grid, theta_grid] = meshgrid(phi, theta); % (nlat x nlon)

% -------------------------------------------------------------------------
%% GRAVITY ANOMALY COMPUTATION (∂V/∂r)
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
deltag_kms2 = GM / R_ref^2 * delta_g;
deltag_mGal = delta_g_kms2 * 1e8;  % 1 km/s^2 = 1e8 mGal

% -------------------------------------------------------------------------
%% CALCULATE BOUGUER ANOMALY
% -------------------------------------------------------------------------
deltag_b = 2*pi*G*rho_crust*elevations; % Bouguer correction
deltag_b_mGal = deltag_b * 1e5; % 1 m/s^2 = 1e5 mGal
scaling_factor = size(deltag_b_mGal,1)/size(deltag_mGal,1);
[a, b] = size(deltag_b_mGal);
if mod(a, scaling_factor) ~= 0 || mod(b, scaling_factor) ~= 0
    error('Matrix dimensions must be divisible by the reduction factor.');
end

% Reshape and average to reduce size of deltag_b
deltag_b_mGal = mean(reshape(deltag_b_mGal, scaling_factor, a/scaling_factor, scaling_factor, b/scaling_factor), [1 3]);
deltag_b_mGal = squeeze(deltag_b_mGal);
BA = deltag_mGal - deltag_b_mGal; % Bouguer Anomaly in mGal

% -------------------------------------------------------------------------
%% PLOTTING
% -------------------------------------------------------------------------
figure;
imagesc(longitudes, latitudes, deltag_b_mGal)
colorbar
figure;
imagesc(longitudes, latitudes, deltag_mGal)
colorbar
figure;
imagesc(longitudes, latitudes, BA)
colorbar