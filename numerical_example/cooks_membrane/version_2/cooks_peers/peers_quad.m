
%% --------------------------------------------------------------------- %%
% ----------------------- MIXED ELEMENT PEERS QUAD ---------------------- %
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
young   = 250 ;                           % modulo di Young
poisson = 4.999 ;                         % modulo di Poisson
ndx     =   4 ;                           % numero suddivisioni in x
ndy     =   2 ;                           % numero suddivisioni in y
nodes   = [0, 0; 48, 44; 48, 60; 0, 44] ;
dl1     = nodes(3,2)-nodes(2,2) ;
dl2     = nodes(4,2) ;
fload   = 100/16;
% Load
f(1,1)  =   0.00 ;                      % load distribiuted direction x
f(2,1)  =   0.00 ;                      % load distribiuted direction y
%
g(1,1) =    0.00 ;                      % traction load direction x edge 1
g(1,2) =    0.00 ;                      % traction load direction y edge 1
%
g(2,1) =    0.00 ;                      % traction load direction x edge 2
g(2,2) =    0.00 ;                      % traction load direction y edge 2
%
g(3,1) =    0.00 ;                      % traction load direction x edge 3
g(3,2) =    0.00 ;                      % traction load direction y edge 3
%
g(4,1) =    0.00 ;                      % traction load direction x edge 4
g(4,2) =   fload ;                      % traction load direction y edge 4
% Boundary conditions (Neumann)
[bn1,bn2,bn3,bn4] = neumann(ndx,ndy,g) ;
bn = [bn1,bn2,bn4] ;
% ----------------------------------------------------------------------- %
lambda = young*poisson/((1+poisson)*(1-2*poisson)) ;
mu = young/(2*(1+poisson)) ;
cf(1,1) = 1/(2*mu) ;
cf(1,2) = -lambda/(4*mu*(mu+lambda)) ;
% ----------------------------------------------------------------------- %

% Geometry
[coordinates,element,mc] = cook(nodes,ndx,ndy,dl1,dl2) ;
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
def = defomesh(spost,element,coordinates) ;

% Plot solution 
plotsol(coordinates,element,stress,spost,rot,def,ndx,ndy) ;
% ----------------------------------------------------------------------- %
