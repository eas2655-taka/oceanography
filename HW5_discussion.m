%% HW5_discussion.m

% safety first
close all
clear all

% include m_map
addpath m_map

% read in TS data
load TS_world_ocean_atlas_2009.mat;

% calculation of dynamic height
g = 9.8; 
p(1,1,:)=depth; % pressure 1dbar = 1m. 

% 3D pressure in dbar
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

%% 1. plotting specific volume at 200 dbar
m=find(depth==200); % find the z-index for 2
figure(1);
m_proj('robinson','clon',180);
m_pcolor(lon,lat,alpha(:,:,m)'); % color shading
hold on;
hold off;
shading flat;
colormap('jet');
caxis([-1 1]*.5e-5);
m_grid('xaxisloc','middle');
m_coast;
colorbar('hori');
title('200 dbar specific volume anomaly based on WOA2009');

%% 2. plotting dynamic height at 0, 200, 600 and 1000 dbar
levs = [0 200 600 1000];

for plt = 1:length(levs);
    
    figure(plt+1);
    m=find(depth==levs(plt));
    m_proj('robinson','clon',180);
    m_pcolor(lon,lat,Z(:,:,m)'); % color shading
    hold on;
    m_contour(lon,lat,Z(:,:,m)',[-1:0.1:3],'k-'); % add contour lines
    hold off;
    shading flat;
    colormap('jet');
    caxis([-1 3]);
    m_grid('xaxisloc','middle');
    m_coast;
    colorbar('hori');
    title([num2str(levs(plt)),' dbar, dynamic height based on WOA2009']);

end

 % Comment: Generally, the magnitude and spatial gradients of dynamic height decreased
 % with higher pressure (deeper depth) levels with the exception of the Southern Ocean. 
 % This means that the geostrophic circulation is generally stronger near the surface. 
