%% plot profile from HOT data

% read data
p=ncread('hotdata.nc','press');
t=ncread('hotdata.nc','theta');
d=ncread('hotdata.nc','mdate');
s=ncread('hotdata.nc','bsal');

% replace missing data (-9) with NaN
t(t==-9)=NaN;
s(s==-9)=NaN;

% plot T-S diagram
x=34:.1:36;
y=-2:1:30;
[x2d,y2d]=meshgrid(x,y);

% calculating density for each point in T-S diagram
rho = sweos(x2d,y2d);
sigma = rho - 1000;

figure(1);
% draw contours
[c,h]=contour(x2d,y2d,sigma);
% contour labeling
clabel(c,h);
ylabel('potential temperature');
xlabel('salinity');

% mmddyy
d(d==-9)=NaN;
mon=floor(d*1e-4);
day=floor((d-mon*1e4)*1e-2);
yr=floor(d-mon*1e4-day*1e2);
dyr = 2000+yr+(mon-1)/12+(day-1)/365.25;

% indices for Dec-Jan-Feb
I=find(mon==12|mon==1|mon==2);

% plot HOT data (S,theta)
hold on;
plot(s(I),t(I),'.');
hold off;

