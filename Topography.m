% Parameters
lmax = 120;
R_ref = 2439.4; % Reference radius in km (from SHADR_HEADER_TABLE)

% Read coefficients
filename = 'gtmes_120v02_sha.tab';
fid = fopen(filename, 'r');
fgetl(fid); % Skip header
data = textscan(fid, '%f%f%f%f%f%f', 'Delimiter', ',', 'CollectOutput', true);
fclose(fid);

coeffs = data{1};
degree = coeffs(:,1);
order  = coeffs(:,2);
Clm    = coeffs(:,3);
Slm    = coeffs(:,4);

% Store in coefficient matrices
Clm_mat = zeros(lmax+1);
Slm_mat = zeros(lmax+1);
for i = 1:length(degree)
    l = degree(i);
    m = order(i);
    Clm_mat(l+1, m+1) = Clm(i);
    Slm_mat(l+1, m+1) = Slm(i);
end

% Grid
nlat = 180;
nlon = 360;
[lon, lat] = meshgrid(linspace(0, 2*pi, nlon), linspace(0, pi, nlat));  % lon = φ, lat = θ
theta = lat(:,1)';  % Colatitude for legendre (1D)

% Initialize topography
r = zeros(nlat, nlon);

% Loop over degrees
for l = 0:lmax
    P_l = legendre(l, cos(theta), 'sch');  % fully normalized (like in PDS description)
    for m = 0:l
        Plm = squeeze(P_l(m+1, :));  % Extract order m
        Plm2D = repmat(Plm', 1, nlon);  % Make into 2D grid

        r = r + Plm2D .* ( ...
            Clm_mat(l+1, m+1) * cos(m * lon) + ...
            Slm_mat(l+1, m+1) * sin(m * lon) );
    end
end

% Final radius and topography
% it should be
% shape_radius(θ, φ) = R_ref + r(θ, φ) from other script
% topography = shape_radius - R_mean

radius = R_ref + r;  % In km
topo = radius - R_ref;  % In km

% Plot
figure;
imagesc(rad2deg(lon(1,:)), 90 - rad2deg(theta), topo); 
axis xy; colorbar;
xlabel('Longitude (deg)');
ylabel('Latitude (deg)');
title('Mercury Topography (m)');


