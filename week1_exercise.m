%% plot profile from HOT data

%% HW1

% read data
p=ncread('hotdata.nc','press');
t=ncread('hotdata.nc','theta');
d=ncread('hotdata.nc','mdate');

% replace missing data (-9) with NaN
t(t==-9)=NaN;

% plot
figure(1);
plot(t,p,'.');
% reverse y axis
set(gca,'Ydir','reverse');
% label axes
xlabel('potential temperature');
ylabel('pressure');

% assemble decimal year
% mmddyy
d(d==-9)=NaN;
mon=floor(d*1e-4);
day=floor((d-mon*1e4)*1e-2);
yr=floor(d-mon*1e4-day*1e2);
dyr = 2000+yr+(mon-1)/12+(day-1)/365.25;

% calculate the surface temperature
sst = t;
sst(p>10)=NaN;

% remove NaNs
dummy = sst+dyr;
sst(isnan(dummy))=[];
dyr(isnan(dummy))=[];

% determine the indices for 2005-2015
I=find(dyr>2005&dyr<2015);

% plot
figure(2);
plot(dyr(I),sst(I),'.-');








