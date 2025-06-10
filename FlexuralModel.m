clear;
close all;
clc;

HOME = pwd;
addpath([HOME '/Data']);
addpath([HOME '/Results']);
addpath([HOME '/Tools'])

%% for Model 3
load([HOME '/Results/Airy_thickness.mat'], 'crust_thickness','longitudes', 'latitudes')
resolution = 1;
[m, n] = size(longitudes);
longitudes = [m, n/64];
[m, n] = size(latitudes);
latitudes = [m, n/64];

%%
latLimT = [-90+(1/resolution/2) 90-(1/resolution/2) 1/resolution]; 
lonLimT = [1/resolution/2  360-(1/resolution/2) 1/resolution]; 

lonT = lonLimT(1):lonLimT(3):lonLimT(2);
latT = fliplr(latLimT(1):latLimT(3):latLimT(2));
LonT = repmat(lonT,length(latT),1);
LatT = repmat(latT',1,length(lonT));

aa = 18;
figure
imagesc(lonT,latT,crust_thickness./1e3);cc=colorbar;
xlabel('Longitude (\circ)','Fontsize',aa)
ylabel('Latitude (\circ)','Fontsize',aa)
ylabel(cc,'Crust Thickness (km)','Fontsize',aa)
set(gca,'YDir','normal','Fontsize',aa)

%% GSHA

cs = GSHA(crust_thickness,179); % it is 1*1 degree plot therefore 179 spherical harmonics
sc = cs2sc(cs);

n = 1:size(sc,1);

D = 100e9*(30e3)^3/(12*(1-0.5^2)); % 200e9 Youngs Modulus, 100e3 Elastic Thickness T_e adapt for Mercury
PHI = (1 + (D)/(500*3.7).*(2.*(n+1)./(2*2439.4e3)).^4).^(-1); %9.81 changed to 3.7, radius changed from 6378 to what is 500?

sc_flex = zeros(size(sc));

for m = 1:size(sc,2)
    sc_flex(:,m) = sc(:,m).*PHI'; %flexural response coeffients adapt to airy function
end

%% GSHS

mapf = GSHS(sc_flex,lonT,latT,179);

%%
figure
imagesc(lonT,latT,mapf./1e3);cc=colorbar;
xlabel('Longitude (\circ)','Fontsize',aa)
ylabel('Latitude (\circ)','Fontsize',aa)
ylabel(cc,'Crust Thickness (km)','Fontsize',aa)
set(gca,'YDir','normal','Fontsize',aa)

figure
imagesc(lonT,latT,(crust_thickness-mapf)./1e3);cc=colorbar;
xlabel('Longitude (\circ)','Fontsize',aa)
ylabel('Latitude (\circ)','Fontsize',aa)
ylabel(cc,'Topography (km)','Fontsize',aa)
set(gca,'YDir','normal','Fontsize',aa)