
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
%
nu = 0.4999;                        % Poisson
mu = 1;                             % Lame constant
lambda = 2*mu*nu/(1-2*nu);          %
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
%[bn1,bn2,bn3,bn4] = neumann(ndx,ndy,g) ;
bn = [] ;
% ----------------------------------------------------------------------- %
cf(1,1) = 1/(2*mu) ;
cf(1,2) = -lambda/(4*mu*(mu+lambda)) ;
% ----------------------------------------------------------------------- %

nl = [2, 4, 8, 16, 32, 64, 128];

name = 'elastic_error_disp_u_peers.txt';
f = fopen(name, 'w');

fprintf(f,'elements v.s. error in L2 norm\n');

for i=1:size(nl,2)

    ndx = nl(i);
    ndy = nl(i);
    % Geometry
    [coordinates,element,mc] = beam(length,heigth,ndx,ndy) ;
    nelem = size(element,1) ; 
    nnod  = size(coordinates,1) ;
    ngdls = max(max(mc)) ;
    ngdd = 2*nelem ;
    ngdr = nnod ;
    ngdlt = ngdls + ngdd + ngdr ;

    % Assembly global system
    [K,load] = assembly_error(coordinates,element,mc,cf,lambda,mu) ; 

    % Solve linear system
    [stress,spost,rot] = solve(K,load,bn,ngdls,ngdd,ngdr) ;

   
    % Compute error in norm L2
    er_u = error_l2_norm(spost, element, coordinates, lambda);    

    % Print results
    fprintf(f, '%2.0f \t %6.5e \n', nelem, er_u);

end
fclose(f);
% ---------------------------------------------------------------------- %
%
% Compute deformate
%def = defomesh(spost,element,coordinates) ;
% Plot solution 
%plotsol(coordinates,element,stress,spost,rot,def,ndx,ndy) ;
 
