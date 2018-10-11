%% HW6_discussion.m

% safety first
close all
clear all

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

%% 1. dynamic height anomaly at 35N
lat0 = 35.5;

% first draw a line at 35.5 N to indicate the location of the transect
figure(1);
m_proj('miller','lon',[280 360],'lat',[0 60]);
m_pcolor(lon,lat,Z(:,:,1)'); % color shading
hold on;
m_contour(lon,lat,Z(:,:,1)',[-1:0.1:3],'k-'); % add contour lines
l=m_line([280 360],[lat0 lat0]);
set(l,'linewidth',3);
set(l,'color','b');
hold off;
shading flat;
colormap('jet');
caxis([-1 3]);
m_grid('xaxisloc','middle');
m_coast;
colorbar('hori');
title('surface dynamic height based on WOA2009');

% then draw the contour map of dynamic height
m=find(lat==lat0); 
figure(2);
pcolor(lon,depth(1:K),squeeze(Z(:,m,:))'); % color shading
set(gca,'ydir','reverse');
hold on;
contour(lon,depth(1:K),squeeze(Z(:,m,:))',-1:.2:3,'k-','linewidth',1);
hold off;
shading flat;
cmp=colormap('jet');
caxis([-1 3]);
axis([280 355 0 1000]);
colorbar('hori');
title(['dynamic height, lat = ',num2str(lat0)]);
xlabel('longitude');
ylabel('deoth');

%% 2. calculate the geostrophic circulation
omega = 2*pi/86400;
f = 2*omega*sin(lat0/180*pi); % coriolis parameter
g = 9.8; % gravity
R = 6.37e6; % earth's radius
dx= (pi/180)*R*cos(pi/180*lat0); % 1degree long in m
dZdx = (Z(2:end,m,:)-Z(1:end-1,m,:))/dx; % calculate dZdx at 30N
lonC = 0.5*(lon(2:end)+lon(1:end-1)); % define long of velocity point
vg = g/f*dZdx; % calculate geostrophic velocity

% plot the meridional (north-south) velocity
figure(3);
pcolor(lonC,depth(1:K),squeeze(vg)'); % color shading
set(gca,'ydir','reverse');
hold on;
contour(lonC,depth(1:K),squeeze(vg)',[0 0],'k-','linewidth',2);
contour(lonC,depth(1:K),squeeze(vg)',[0.01:0.01:0.1],'k-','linewidth',1);
contour(lonC,depth(1:K),squeeze(vg)',[-0.1:0.01:-0.01],'k--','linewidth',1);
hold off;
shading flat;
cmp=colormap('jet');
caxis([-0.1 0.1]);
axis([280 355 0 1000]);
colorbar('hori');
title(['northward geostrophic velocity, lat = ',num2str(lat0)]);
xlabel('longitude');
ylabel('deoth');

