
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
%                          b( v, div(Sigma) )               = -( f,v )    %
%                          c( psi, as(sigma) )              = 0           %
% Sigma, Tau are tensor   (2x2)                                           %
% u, v       are vector   (2x1)                                           %
% eta, psi   are costant  (1x1)                                           %
% ----------------------------------------------------------------------- %

%% INPUT DATI 
clear; close all; clc;
% Geometry
length  =  1 ;                      % lunghezza trave
heigth  =  1 ;                      % altezza trave
%
%poisson = 0.4999;                   % modulo di Poisson
%G = 1;                              % Lame constant mu = Shear modulus   

G = 79.3;
lambda = 123;

%
g(1,1) =   0.00 ;                   % traction load direction x edge 1
g(1,2) =   0.00 ;                   % traction load direction y edge 1
%
g(2,1) =   0.00 ;                   % traction load direction x edge 2
g(2,2) =   0.00 ;                   % traction load direction y edge 2
%
g(3,1) =   0.00 ;                   % traction load direction x edge 3
g(3,2) =   0.00 ;                   % traction load direction y edge 3
%
g(4,1) =   0.00 ;                   % traction load direction x edge 4
g(4,2) =   0.00 ;                   % traction load direction y edge 4
% Boundary conditions (Neumann)
bn = [] ;
% ----------------------------------------------------------------------- %
%young  =  G*(1+poisson)*2;                    % modulo di Young
%lambda = young*poisson/((1+poisson)*(1-2*poisson)) ;
% ----------------------------------------------------------------------- %

nl = [2, 4, 8, 16, 32, 64, 128];

name_u = 'elastic_error_disp_u_peers_2bm.txt';
name_r = 'elastic_error_rot_u_peers_2bm.txt';

f  = fopen(name_u, 'w');
ff = fopen(name_r, 'w');

fprintf(f,'elements v.s. error in L2 norm\n');
fprintf(ff,'elements v.s. error in L2 norm\n');


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
    [K,load] = assembly_error(coordinates,element,mc,lambda,G) ; 

    % Solve linear system
    [stress,spost,rot] = solve(K,load,bn,g,ndx,ndy,ngdls,ngdd,ngdr) ;

    % Compute error in norm L2
    er_u = error_l2_norm(spost, mc, element, coordinates, lambda);    

    % Compute error rotations in norm L2
    er_r = error_l2_norm_rot(rot, mc, element, coordinates, lambda);
    
    % Print results
    fprintf(f, '%2.0f \t %6.5e \n', nelem, er_u);
    fprintf(ff, '%2.0f \t %6.5e \n', nelem, er_r);

end
fclose(f);
fclose(ff);
% Compute deformate
%def = defomesh(spost,element,coordinates) ;
% Plot solution 
%plotsol(coordinates,element,stress,spost,rot,def,ndx,ndy) ;
% ---------------------------------------------------------------------- %
