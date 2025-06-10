% run_gshs_example.m
clear;
close all;
clc;

HOME = pwd;
addpath([HOME '/Data']);
addpath([HOME '/Results']);
addpath([HOME '/Tools'])

% Load the spherical harmonics coefficients from the file
filename = [HOME '/Data/ggmes_50v06_sha.tab'];
maxDegree = 100;
R_ref = 2439.4e3;                     % Reference radius (m)
GM = 22031.815e9;              % Mercury GM (m^3/s^2)
coeffs = readmatrix(filename, 'FileType', 'text', 'Delimiter', ',');
resolution = 1;

% Determine maximum degree (lmax)
lmax = max(coeffs(:,1));  % Assuming column 1 = degree n

% Initialize field matrix (sc format)
field = zeros(lmax+1, 2*lmax+1);

% Fill the field matrix from the data
for i = 1:size(coeffs,1)
    n = coeffs(i,1);  % degree
    m = coeffs(i,2);  % order
    C = coeffs(i,3);  % Cnm
    S = coeffs(i,4);  % Snm
    
    field(n+1, lmax+1+m) = C;  % column index for Cnm = lmax+1 + m
    if m ~= 0
        field(n+1, lmax+1-m) = S;  % column index for Snm = lmax+1 - m
    end
end

% Define the latitude and longitude grid
latLimT = [-90+(1/resolution/2) 90-(1/resolution/2) 1/resolution]; 
lonLimT = [1/resolution/2  360-(1/resolution/2) 1/resolution]; 

lonT = lonLimT(1):lonLimT(3):lonLimT(2);
latT = fliplr(latLimT(1):latLimT(3):latLimT(2));
%latT = latLimT(1):latLimT(3):latLimT(2);
LonT = repmat(lonT,length(latT),1);
LatT = repmat(latT',1,length(lonT));

% Run the spherical harmonic synthesis
f = GSHS(field, lonT, 90-latT, maxDegree);  % Output is [length(th) x length(lon)]


% Plot the result
figure
aa = 18;
imagesc(lonT,latT,delta_g_mGal);cc=colorbar;
xlabel('Longitude (\circ)','Fontsize',aa)
ylabel('Latitude (\circ)','Fontsize',aa)
ylabel(cc,'Gravity Anomaly (mGal)','Fontsize',aa)
set(gca,'YDir','normal','Fontsize',aa)
