function density = swdens(S,T,P)

%% input: S, T, P
%% S = salinity, psu
%% T = temperature, deg C
%% P = pressure, db
%% output: density of seawater, kg/m3, based on UNESCO 1981 equation of state

%% set up coefficient
coeff = [9.99843E+02 6.79395E-02 -9.09529E-03 1.00169E-04 -1.12008E-06 6.53633E-09; 
         8.24493E-01 -4.08990E-03 7.64380E-05 -8.24670E-07 5.38750E-09 0.00000E+00;
        -5.72466E-03 1.02270E-04 -1.65460E-06 0.00000E+00 0.00000E+00 0.00000E+00;
        4.83140E-04 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00	0.00000E+00];

%% calculate freshwater density
dens0 = coeff(1,1)+coeff(1,2)*T+coeff(1,3)*T.^2+... 
    coeff(1,4)*T.^3+coeff(1,5)*T.^4+coeff(1,6)*T.^5;

%% salinity correction
Scorr = S.*(coeff(2,1)+coeff(2,2)*T+coeff(2,3)*T.^2+... 
    coeff(2,4)*T.^3+coeff(2,5)*T.^4+coeff(2,6)*T.^5)+...
    S.^(3/2).*(coeff(3,1)+coeff(3,2)*T+coeff(3,3)*T.^2+... 
    coeff(3,4)*T.^3+coeff(3,5)*T.^4+coeff(3,6)*T.^5)+...
    S.^2*coeff(4,1);

%% density at surface pressure
dens0 = dens0 + Scorr;

% COMPUTE COMPRESSION TERMS
P = P/10;  %convert from db to atmospheric pressure units
T68 = T * 1.00024;

% Pure water terms of the secant bulk modulus at atmos pressure.
% UNESCO eqn 19 p 18
h3 = -5.77905E-7;
h2 = +1.16092E-4;
h1 = +1.43713E-3;
h0 = +3.239908;   %[-0.1194975];
AW  = h0 + (h1 + (h2 + h3.*T68).*T68).*T68;
k2 =  5.2787E-8;
k1 = -6.12293E-6;
k0 =  +8.50935E-5;   %[+3.47718E-5];
BW  = k0 + (k1 + k2*T68).*T68;
e4 = -5.155288E-5;
e3 = +1.360477E-2;
e2 = -2.327105;
e1 = +148.4206;
e0 = 19652.21;    %[-1930.06];
KW  = e0 + (e1 + (e2 + (e3 + e4*T68).*T68).*T68).*T68;   % eqn 19

% SEA WATER TERMS OF SECANT BULK MODULUS AT ATMOS PRESSURE.
j0 = 1.91075E-4;
i2 = -1.6078E-6;
i1 = -1.0981E-5;
i0 =  2.2838E-3;
SR = sqrt(S);
A  = AW + (i0 + (i1 + i2*T68).*T68 + j0*SR).*S;
m2 =  9.1697E-10;
m1 = +2.0816E-8;
m0 = -9.9348E-7;
B = BW + (m0 + (m1 + m2*T68).*T68).*S;   % eqn 18
f3 =  -6.1670E-5;
f2 =  +1.09987E-2;
f1 =  -0.603459;
f0 = +54.6746;
g2 = -5.3009E-4;
g1 = +1.6483E-2;
g0 = +7.944E-2;
K0 = KW + (  f0 + (f1 + (f2 + f3*T68).*T68).*T68 ...
        +   (g0 + (g1 + g2*T68).*T68).*SR         ).*S;      % eqn 16
K = K0 + (A + B.*P).*P;  % eqn 15

%% density including the compression effect
density = dens0./(1-P./K);
return
