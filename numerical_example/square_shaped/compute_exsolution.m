% Compute load starting to the solutions
clear; clc;
syms x y lm mu

u = -2.*x.^2.*y.*(1-x).^2.*(1-3.*y+2.*y.^2);
v =  2.*x.*y.^2.*(1-y).^2.*(1-3.*x+2.*x.^2);
%u = cos(pi.*x).*sin(2.*pi.*y);
%v = sin(pi.*x).*cos(pi.*y);

e11 = diff(u,x);
e22 = diff(v,y);
e12 = (diff(u,y)+diff(v,x)).*0.5;
e21 = e12;


sig11 = 2.*mu.*e11 + lm.*(e11 + e22);
sig22 = 2.*mu.*e22 + lm.*(e11 + e22);
sig12 = 2.*mu.*e12;
sig21 = sig12;


rot =  simplify( (diff(u,y)-diff(v,x))/2 );
f1 = simplify( diff(sig11,x) + diff(sig12,y) );
f2 = simplify( diff(sig21,x) + diff(sig22,y) );
