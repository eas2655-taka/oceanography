%% week5_exercise.m

% safety first
close all
clear all

% include m_map
addpath m_map

% read in TS data
load TS_world_ocean_atlas_2009.mat;

% practice making a map
% figure(1);
% m_proj('robinson','clon',180);
% m_pcolor(lon,lat,temperature(:,:,1)');
% shading flat;
% m_grid;
% m_coast;

% calculation of dynamic height
g = 9.8; 
p(1,1,:)=depth;
% 3D pressure
p3d=repmat(p,[360 180 1]);
% calculating density
rho=swdens(salinity,temperature,p3d);
% calculating ref density
rhoref=swdens(35,0,p3d);
% calculating alpha (specific volume anomaly)
alpha=1./rho - 1./rhoref;

% now vertically integrate alpha to get Z (dynamic topo)

% ref pressure = 2000 dbar
K = find(depth==2000);

% use for loop to integrate
Z=zeros(360,180,K);
for k = K:-1:2
    Z(:,:,k-1)=Z(:,:,k)+alpha(:,:,k)/g*(depth(k)-depth(k-1))*1e4;
end

% plotting dynamic height
figure(2);
m_proj('robinson','clon',180);
m_pcolor(lon,lat,Z(:,:,1)'); % color shading
hold on;
m_contour(lon,lat,Z(:,:,1)',[-1:0.1:3],'k-'); % add contour lines
hold off;
shading flat;
colormap('jet');
caxis([-1 3]);
m_grid;
m_coast;
colorbar('hori');
title('surface dynamic height based on WOA2009');

% plotting specific volume
figure(3);
m_proj('robinson','clon',180);
m_pcolor(lon,lat,alpha(:,:,1)'); % color shading
hold on;
%m_contour(lon,lat,Z(:,:,1)',[-1:0.1:3],'k-'); % add contour lines
hold off;
shading flat;
colormap('jet');
caxis([-1 1]*1e-5);
m_grid;
m_coast;
colorbar('hori');
title('surface specific volume anomaly based on WOA2009');

% plotting in polar projection
figure(4);
m_proj('stereographic','lat',-90,'radius',60); % polar projection
m_pcolor(lon,lat,Z(:,:,1)'); % color shading
hold on;
m_contour(lon,lat,Z(:,:,1)',[-1:0.1:3],'k-'); % add contour lines
hold off;
shading flat;
colormap('jet');
caxis([-1 3]);
m_grid('xaxisloc','top');
m_coast;
colorbar('hori');
title('surface dynamic height based on WOA2009');










