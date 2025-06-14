clear;
close all;
clc;

HOME = pwd;
addpath([HOME '/Data']);
addpath([HOME '/Results']);
addpath([HOME '/Tools'])
addpath([HOME '/ScientificColourMaps8/vik']);
load('vik.mat');

%% for Model 3
load([HOME '/Results/Airy_thickness.mat'], 'crust_thickness');
load([HOME '/Results/coeffs_obs.mat'], 'V');
load([HOME '/Results/elevations.mat'], 'elevations');
resolution = 1;
Te = 30e3;

%%
latLim = [-90+(1/resolution/2) 90-(1/resolution/2) 1/resolution]; 
lonLim = [(1/resolution/2) 360-(1/resolution/2) 1/resolution]; 

lonT = lonLim(1):lonLim(3):lonLim(2);
latT = latLim(1):latLim(3):latLim(2);
LonT = repmat(lonT,length(latT),1);
LatT = repmat(latT',1,length(lonT));

aa = 16;
% figure
% imagesc(lonT,latT,crust_thickness./1e3);cc=colorbar;
% xlabel('Longitude (\circ)','Fontsize',aa)
% ylabel('Latitude (\circ)','Fontsize',aa)
% ylabel(cc,'Airy Crust Thickness (km)','Fontsize',aa)
% set(gca,'YDir','normal','Fontsize',aa)

%% GSHA

cs = GSHA(crust_thickness,50); % it is 1*1 degree plot therefore 179 spherical harmonics
sc = cs2sc(cs);

n = 1:size(sc,1);
nu = 0.25;
R = 2439.4e3;
D = 100e9*(Te)^3/(12*(1-nu^2)); % 200e9 Youngs Modulus, 100e3 Elastic Thickness T_e and oissons Ratio adapt for Mercury
%PHI = (1 + (D)/(400*3.7).*(2.*(n+1)./(2*2439.4e3)).^4).^(-1); %9.81 changed to 3.7, radius changed from 6378 to what is 500 density difference 
PHI = (1+D/(400*3.7).*(1/R.^4.*(n.*(n+1)-2).^2./(1-(1-nu.^2)./(n.*(n+1)))+12.*(1-nu)./(Te.^2.*R.^2)*(1-2./(n.*(n+1))./(1-(1-nu)./(n.*(n+1)))))).^(-1);
% between mantle and crust

sc_flex = zeros(size(sc));

for m = 1:size(sc,2)
    sc_flex(:,m) = sc(:,m).*PHI'; %flexural response coeffients adapt to airy function
end

%% GSHS

mapf = GSHS(sc_flex,lonT,90-latT,50);

%%
figure
imagesc(lonT,latT,mapf./1e3);cc=colorbar;
set(gca, 'YDir', 'normal','Fontsize', aa);          % Flip Y so North is up
xlabel('Longitude ($^\circ$)', 'Interpreter', 'latex', 'Fontsize', aa)
ylabel('Latitude ($^\circ$)', 'Interpreter', 'latex', 'Fontsize', aa)
c = colorbar;
ylabel(c, 'Flexural Crust Thickness (km)', 'Interpreter', 'latex', 'Fontsize', aa);
c = colormap(vik);
set(gca, 'ylim', [-90 90]);
set(gca, 'ytick', -90:30:90);
set(gca, 'xlim', [0 360]);
set(gca, 'xtick', 0:30:360);
axis equal tight;
saveas(gcf, 'Figures/Flexural_Thickness.svg');


figure
imagesc(lonT,latT,(crust_thickness-mapf)./1e3);cc=colorbar;


%% save
% save also V_Model and insert it directly
save([HOME '/Results/Flexural_thickness.mat'], 'mapf', 'lonT', 'latT');
