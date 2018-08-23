% plot profile from HOT data

%% HW1

% read data
p=ncread('hotdata.nc','press');
t=ncread('hotdata.nc','theta');

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

