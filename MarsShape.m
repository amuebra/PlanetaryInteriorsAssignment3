clear;
close all;
clc;

filename = 'mega90n000cb.img';
resolution = 4;

% Read in the file.
f = fopen(filename,'r','ieee-be');
el1 = fread(f,[360*resolution Inf],'int16')';
fclose(f); 

filename = 'megc90n000cb.img';
resolution = 4;

% Read in the file.
f = fopen(filename,'r','ieee-be');
el2 = fread(f,[360*resolution Inf],'int16')';
fclose(f); 

filename = 'megr90n000cb.img';
resolution = 4;

% Read in the file.
f = fopen(filename,'r','ieee-be');
el3 = fread(f,[360*resolution Inf],'int16')';
fclose(f); 

filename = 'megt90n000cb.img';
resolution = 4;

% Read in the file.
f = fopen(filename,'r','ieee-be');
el4 = fread(f,[360*resolution Inf],'int16')';
fclose(f); 

%%

latLimT = [-90+(1/resolution/2) 90-(1/resolution/2) 1/resolution]; 
lonLimT = [1/resolution/2  360-(1/resolution/2) 1/resolution]; 

lonT = lonLimT(1):lonLimT(3):lonLimT(2);
latT = fliplr(latLimT(1):latLimT(3):latLimT(2));
LonT = repmat(lonT,length(latT),1);
LatT = repmat(latT',1,length(lonT));

load('vik.mat')
aa = 18;

figure
subplot(2,2,1)
imagesc(lonT,latT,el1./1e3);cc=colorbar;colormap(vik);
xlabel('Longitude (\circ)','Fontsize',aa)
ylabel('Latitude (\circ)','Fontsize',aa)
ylabel(cc,'areoide (km)','Fontsize',aa)
set(gca,'YDir','normal','Fontsize',aa)
subplot(2,2,2)
imagesc(lonT,latT,el2);cc=colorbar;colormap(vik);
xlabel('Longitude (\circ)','Fontsize',aa)
ylabel('Latitude (\circ)','Fontsize',aa)
ylabel(cc,'count (-)','Fontsize',aa)
set(gca,'YDir','normal','Fontsize',aa)
subplot(2,2,3)
imagesc(lonT,latT,el3./1e3);cc=colorbar;colormap(vik);
xlabel('Longitude (\circ)','Fontsize',aa)
ylabel('Latitude (\circ)','Fontsize',aa)
ylabel(cc,'radius offset to 3396 (km)','Fontsize',aa)
set(gca,'YDir','normal','Fontsize',aa)
subplot(2,2,4)
imagesc(lonT,latT,el4./1e3);cc=colorbar;colormap(vik);
xlabel('Longitude (\circ)','Fontsize',aa)
ylabel('Latitude (\circ)','Fontsize',aa)
ylabel(cc,'topography (km)','Fontsize',aa)
set(gca,'YDir','normal','Fontsize',aa)

