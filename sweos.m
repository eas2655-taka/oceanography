function density = sweos(S,T)

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

%% total density
density = dens0 + Scorr;
