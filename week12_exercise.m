%% 2 box model of carbon

close all
clear all

% parameters
A = 4e14;   % m2, surf area
Vs = A*100; % m3, surf box vol
Vd = A*3000; % m3, deep box
K0 = 0.04;   % solubility of CO2
K1 = 1e-6;   % K1
K2 = 1e-9;   % K2
Ps = 1.0e-6; % mol/L, surface P
Pd = 2.2e-6; % mol/L, deep P
Vm = 3e8;    % m3/s
G  = 2000/3e7; % m/s, gas transfer coefficient
Rcp= 106;    %C:P ratio

% Initial condition
Cs(1)=2000e-6; % mol/L, surface C
Cd(1)=2200e-6; % mol/L, surface C
pCO2atm = 280e-6; % atm, atmos CO2
Alk = 2300e-6; % mol/L, surf alkalinity

% time stepping parameters
dt = 60*60*24; % 1 day in sec
N  = 365*1000; % 1000 year integration
year(1)=0;

% time stepping loop
for i=1:N-1 
    % calculate ocean pCO2
    pCO2ocn(i)=K2*(2*Cs(i)-Alk)^2/(K0*K1*(Alk-Cs(i)));
    
    % transport terms
    Gasex(i)=-G*A*K0*(pCO2ocn(i)-pCO2atm);
    Circ(i)=Vm*(Cd(i)-Cs(i));
    Bio(i)=Vm*(Pd-Ps)*Rcp;
    
    % rate of change
    dCsdt=1/Vs*(Gasex(i)+Circ(i)-Bio(i));
    dCddt=1/Vd*(-Circ(i)+Bio(i));
    
    % step forward in time
    Cs(i+1)=Cs(i)+dt*dCsdt;
    Cd(i+1)=Cd(i)+dt*dCddt;
    
    % time
    year(i+1)=year(i)+dt/(60*60*24*365);
end

% diagnostics
HCO3 = 2*Cs - Alk;
CO3  = Alk-Cs;
H = K2*HCO3./CO3;
pH=-log10(H);

% plot output
figure(1);
subplot(2,1,1);
plot(year,Cs*1e6,'b-');
hold on;
plot(year,Cd*1e6,'k-');
hold off;
ylabel('DIC, micro-M');
xlabel('time, year');

subplot(2,1,2);
plot(year,pH,'b-');
ylabel('pH');
xlabel('time, year');




