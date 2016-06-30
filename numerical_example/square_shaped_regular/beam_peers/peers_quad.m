
%% --------------------------------------------------------------------- %%
% ----------------------- MIXED ELEMENT PEERS QUAD ---------------------- %
%                      FOR LINEAR ELASTICITY PROBLEM                      %
% ---------------   ( by Marco Pingaro & Daniele Boffi )   -------------- %
%% Phd Student                                                             %
% mail  : marco.pingaro@iusspavia.it                                      %
% ----------------------------------------------------------------------- %
% The problem is : laplace(u) = f                                         %
% Dirichlet condition is u_D = 0 on boundary domain                       %
%                                                                         %
% The Weak Formulation is:                                                %
% a(C^(-1)Sigma,Tau ) + b( u,div(Tau) ) + c( eta, as(tau) ) = 0           %
%                          b( v, div(Sigma) )               = ( f,v )     %
%                          c( psi, as(sigma) )              = 0           %
% Sigma, Tau are tensor   (2x2)                                           %
% u, v       are vector   (2x1)                                           %
% eta, psi   are costant  (1x1)                                           %
% ----------------------------------------------------------------------- %

%% INPUT DATI 
clear all; close all; clc;
% Geometry
length  =   1;                      % lunghezza trave
heigth  =   1 ;                     % altezza trave

young   = 1.0 ;                     % modulo di Young
poisson = 0.0 ;                     % modulo di Poisson
%lambda = 123.0;                     % Lame constant
%mu = 79.3;                          % Lame constant

ndx     =   1 ;                     % numero suddivisioni in x
ndy     =   1 ;                     % numero suddivisioni in y
% Load
f(1,1)  =  1.0 ;                    % load distribiuted direction x
f(2,1)  =  1.0 ;                    % load distribiuted direction y
%
g(1,1) =   0.0 ;                    % traction load direction x edge 1
g(1,2) =   0.0 ;                    % traction load direction y edge 1
%
g(2,1) =   0.0 ;                    % traction load direction x edge 2
g(2,2) =   0.0 ;                    % traction load direction y edge 2
%
g(3,1) =   0.0 ;                    % traction load direction x edge 3
g(3,2) =   0.0 ;                    % traction load direction y edge 3
%
g(4,1) =   0.0 ;                    % traction load direction x edge 4
g(4,2) =   0.0 ;                    % traction load direction y edge 4
% Boundary conditions (Neumann)
[bn1,bn2,bn3,bn4] = neumann(ndx,ndy,g) ;
bn = [] ;
% ----------------------------------------------------------------------- %


lambda = young*poisson/((1+poisson)*(1-2*poisson)) ;
mu = young/(2*(1+poisson)) ;
cf(1,1) = 1/(2*mu) ;
cf(1,2) = -lambda/(4*mu*(mu+lambda)) ;
% ----------------------------------------------------------------------- %

s = [2, 4, 8, 16, 32, 64, 128];

name_1 = 'elastic_error_disp_ux.txt';
name_2 = 'elastic_error_disp_uy.txt';
f1 = fopen(name_1, 'w');
f2 = fopen(name_2, 'w');

fprintf(f1,'elements v.s. error in L2 norm\n');
fprintf(f2,'elements v.s. error in L2 norm\n');

for i=1:3

    ndx = s(i);
    ndy = s(i);
    % Geometry
    [coordinates,element,mc] = beam(length,heigth,ndx,ndy) ;
    nelem = size(element,1) ; 
    nnod  = size(coordinates,1) ;
    ngdls = max(max(mc)) ;
    ngdd = 2*nelem ;
    ngdr = nnod ;
    ngdlt = ngdls + ngdd + ngdr ;

    % Assembly global system
    [K,load] = assembly(coordinates,element,mc,cf,f) ; 

    % Solve linear system
    [stress,spost,rot] = solve(K,load,bn,ngdls,ngdd,ngdr) ;

    % Compute deformate
    %def = defomesh(spost,element,coordinates) ;

    % Plot solution 
    %plotsol(coordinates,element,stress,spost,rot,def,ndx,ndy) ;
    
    % Compute error in norm L2
    [er_ux, er_uy] = error_l2_norm(spost, element, coordinates);    

    % Print results
    fprintf(f1, '%6.4f \t %6.5e \n', ndx, er_ux);
    fprintf(f2, '%6.4f \t %6.5e \n', ndx, er_uy);

end
fclose(f1);
fclose(f2);
% ----------------------------------------------------------------------- %
