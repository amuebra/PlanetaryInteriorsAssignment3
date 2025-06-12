clear;
close all;
clc;

HOME = pwd;
addpath([HOME '/Data']);
addpath([HOME '/Results']);
addpath([HOME '/Tools'])

%% for Model 3
load([HOME '/Results/Airy_thickness.mat'], 'crust_thickness')
resolution = 1;

%%
latLimT = [-90+(1/resolution/2) 90-(1/resolution/2) 1/resolution]; 
lonLimT = [(1/resolution/2) 360-(1/resolution/2) 1/resolution]; 

lonT = lonLimT(1):lonLimT(3):lonLimT(2);
latT = latLimT(1):latLimT(3):latLimT(2);
LonT = repmat(lonT,length(latT),1);
LatT = repmat(latT',1,length(lonT));

aa = 18;
figure
imagesc(lonT,latT,crust_thickness./1e3);cc=colorbar;
xlabel('Longitude (\circ)','Fontsize',aa)
ylabel('Latitude (\circ)','Fontsize',aa)
ylabel(cc,'Airy Crust Thickness (km)','Fontsize',aa)
set(gca,'YDir','normal','Fontsize',aa)

%% GSHA

cs = GSHA(crust_thickness,179); % it is 1*1 degree plot therefore 179 spherical harmonics
sc = cs2sc(cs);

n = 1:size(sc,1);

D = 100e9*(30e3)^3/(12*(1-0.25^2)); % 200e9 Youngs Modulus, 100e3 Elastic Thickness T_e and oissons Ratio adapt for Mercury
PHI = (1 + (D)/(400*3.7).*(2.*(n+1)./(2*2439.4e3)).^4).^(-1); %9.81 changed to 3.7, radius changed from 6378 to what is 500 density difference 
% between mantle and crust

sc_flex = zeros(size(sc));

for m = 1:size(sc,2)
    sc_flex(:,m) = sc(:,m).*PHI'; %flexural response coeffients adapt to airy function
end

%% GSHS

mapf = GSHS(sc_flex,lonT,90-latT,179);

%%
figure
imagesc(lonT,latT,mapf./1e3);cc=colorbar;
xlabel('Longitude (\circ)','Interpreter', 'Latex', 'Fontsize',aa)
ylabel('Latitude (\circ)','Interpreter', 'Latex','Fontsize',aa)
ylabel(cc,'Flexural Crust Thickness (km)','Interpreter', 'Latex','Fontsize',aa)
set(gca,'YDir','normal','Fontsize',aa)

figure
imagesc(lonT,latT,(crust_thickness-mapf)./1e3);cc=colorbar;
xlabel('Longitude (\circ)','Interpreter', 'Latex','Fontsize',aa)
ylabel('Latitude (\circ)','Interpreter', 'Latex','Fontsize',aa)
ylabel(cc,'Residual between Airy and Flexural (km)','Interpreter', 'Latex','Fontsize',aa)
set(gca,'YDir','normal','Fontsize',aa)

%% save
save([HOME '/Results/Flexural_thickness.mat'], 'mapf', 'lonT', 'latT');