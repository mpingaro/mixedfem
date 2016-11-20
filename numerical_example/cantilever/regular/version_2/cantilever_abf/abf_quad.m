
%% --------------------------------------------------------------------- %%
% --------------------- MIXED ELEMENT PEERS-ABF QUAD -------------------- %
%                      FOR LINEAR ELASTICITY PROBLEM                      %
% ---------------   ( by Marco Pingaro & Daniele Boffi )   -------------- %
% Phd Student                                                             %
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
length  = 10 ;                     % lunghezza trave
heigth  =  2 ;                     % altezza trave
young   = 1500 ;                   % modulo di Young
poisson = 0.4999 ;                 % modulo di Poisson
% 
g(1,1) = 0.0 ;                     % traction load direction x edge 1  
g(1,2) = 0.0 ;                     % traction load direction y edge 1
%
g(2,1) = 0.0 ;                     % traction load direction x edge 2
g(2,2) = 0.0 ;                     % traction load direction y edge 2
%
g(3,1) = 0.0 ;                     % traction load direction x edge 3
g(3,2) = 0.0 ;                     % traction load direction y edge 3
%
g(4,1) = 300.0 ;                   % traction load direction x edge 4
g(4,2) = 0.0 ;                     % traction load direction y edge 4


nx = [4,8,16,32,64,128];
ny = [2,4, 8,16,32, 64];

fname = 'error_beam_u_l2_abf_ver2.txt';
f = fopen( fname, 'w');
fprintf(f, 'element vs. error u in norm L2\n');

for i=1:numel(nx)

ndx = nx(i);
ndy = ny(i);

% ----------------------------------------------------------------------- %
lambda = young*poisson/((1+poisson)*(1-2*poisson)) ;
mu = young/(2*(1+poisson)) ;
cf(1,1) = 1/(2*mu) ;
cf(1,2) = -lambda/(4*mu*(mu+lambda)) ;

% Geometry
[coordinates,element,mc] = beam(length,heigth,ndx,ndy) ;
nelem = size(element,1) ; 
nnod  = size(coordinates,1) ;
ngdls = max(max(mc)) ;
ngdd = 6*nelem ;
ngdr = nnod ;
ngdlt = ngdls + ngdd + ngdr ;

% Boundary conditions
[bn1,bn2,bn3,bn4] = neumann(coordinates,ndx,ndy,g,heigth) ;
bn = [bn1,bn2,bn3,bn4] ;

% Assembly global system
[K,load] = assembly(coordinates,element,mc,cf) ; 

% Solve linear system
[stress,spost,rot] = solve(K,load,bn,ngdls,ngdd,ngdr) ;

% Compute deformate
%def = defomesh(spost,element,coordinates) ;

% Plot solution
%plotsol(coordinates,element,stress,spost,rot,def,ndx,ndy) ;

%
er_u = error_beam_l2_norm(spost,element,coordinates,lambda,young,poisson,g(4,1));
fprintf(f, '%6.0f \t %6.5e \n', nelem, er_u);

end
fclose(f);
% ------------------------------------------------------------------------ %
