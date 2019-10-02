function [fnc_fin,fnc_high,fnc_low,lAl] = fit_condest_rho(A,Q)
%FIT_CONDEST_RHO Summary of this function goes here
%   Detailed explanation goes here
lin_coeff=@(x1,x2,y1,y2)[(y1 - y2)/(x1 - x2); (x1*y2 - x2*y1)/(x1 - x2)];

lAl = eigs(A, 1);
rho0=lAl*1000;
rho1=lAl*10000;
rho2=lAl/10;
rho3=lAl/100;

y0=log10(eigs(A+rho0*Q, 1));
x0=log10(rho0);
y1 = log10(eigs(A+rho1*Q, 1));
x1=log10(rho1);
coeff0=lin_coeff(x0,x1,y0,y1);
fnc_high=@(x) coeff0(1)*x+coeff0(2);

y2 = log10(eigs(A+rho2*Q, 1));
y3 = log10(eigs(A+rho3*Q, 1));
x2 = log10(rho2);
x3=log10(rho3);


coeff=lin_coeff(x2,x3,y2-fnc_high(x2)-log10(lAl),y3-fnc_high(x3),-log10(lAl));

fnc_low=@(x) coeff(1)*x+coeff(2);
fnc_fin=@(x)fnc_high(x)+fnc_low(x);

end

