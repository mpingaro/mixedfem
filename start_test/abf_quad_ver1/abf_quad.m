
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
length  = 4 ;                      % lunghezza trave
heigth  = 1 ;                      % altezza trave
young   = 1.0 ;                    % modulo di Young
poisson = 0.3 ;                    % modulo di Poisson
ndx     =  4 ;                    % numero suddivisioni in x
ndy     =  1 ;                    % numero suddivisioni in y
% Load
f(1,1) =  0.00 ;                   % load distribiuted direction x
f(2,1) =  1.00 ;                   % load distribiuted direction y
% 
g(1,1) =  0.00 ;                   % traction load direction x edge 1  
g(1,2) =  0.00 ;                   % traction load direction y edge 1
%
g(2,1) =  0.00 ;                   % traction load direction x edge 2
g(2,2) =  0.00 ;                   % traction load direction y edge 2
%
g(3,1) =  0.00 ;                   % traction load direction x edge 3
g(3,2) =  0.00 ;                   % traction load direction y edge 3
%
g(4,1) =  0.00 ;                   % traction load direction x edge 4
g(4,2) =  0.00 ;                   % traction load direction y edge 4
% Boundary conditions
bn = 3 ;
% ----------------------------------------------------------------------- %
lambda = young*poisson/((1+poisson)*(1-2*poisson)) ;
G = young/(2*(1+poisson)) ;
% Geometry
[coordinates,element,mc] = beam(length,heigth,ndx,ndy) ;
nelem = size(element,1) ; 
nnod  = size(coordinates,1) ;
ngdls = max(max(mc)) ;
ngdd = 6*nelem ;
ngdr = nnod ;
ngdlt = ngdls + ngdd + ngdr ;

% Assembly global system
[K,load] = assembly(coordinates,element,mc,lambda,G,f) ;  

% Solve linear system
[stress,spost,rot] = solve(K,load,bn,g,ndx,ndy,ngdls,ngdd,ngdr) ;

% Compute deformate
def = defomesh(spost,element,coordinates) ;

% Plot solution
%plotsol(coordinates,element,stress,spost,rot,def,ndx,ndy) ;
% ------------------------------------------------------------------------ %
