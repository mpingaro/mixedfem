
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
clear; close all; clc;
% Geometry
length  = 1 ;                      % lunghezza trave
heigth  = 1 ;                      % altezza trave
%nu = 0.4999;                       % Poisson
%G  = 1.;                           % Lame constant = Shear modulus
%lambda = 2*G*nu/(1-2*nu);         % Lame constant
G = 79.3;
lambda = 123;
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
bn = [] ;
% ----------------------------------------------------------------------- %
nl = [2, 4, 8, 16, 32, 64];

name_u = 'elastic_error_disp_u_abf_trapezioidal_1b.txt';
name_r = 'elastic_error_rot_u_abf_trapezoidal_1b.txt';
name_div = 'elastic_error_div_u_abf_trapezoidal_1b.txt';

f  = fopen(name_u, 'w');
ff = fopen(name_r, 'w');
fff= fopen(name_div, 'w');

fprintf(f,'elements v.s. error in L2 norm\n');
fprintf(ff,'elements v.s. error in L2 norm\n');
fprintf(fff,'elements v.s. error in L2 norm\n');

for i=1:size(nl,2)
    
    ndx = nl(i);
    ndy = nl(i);

    % Geometry
    [coordinates,element,mc] = beam(length,heigth,ndx,ndy) ;
    nelem = size(element,1) ; 
    nnod  = size(coordinates,1) ;
    ngdls = max(max(mc)) ;
    ngdd = 6*nelem ;
    ngdr = nnod ;
    ngdlt = ngdls + ngdd + ngdr ;

    % Assembly global system
    [K,load] = assembly_error(coordinates,element,mc,lambda,G) ;  

    % Solve linear system
    [stress,spost,rot] = solve(K,load,bn,g,ndx,ndy,ngdls,ngdd,ngdr) ;

    % Compute error in norm L2
    er_u = error_l2_norm(spost, element, coordinates, lambda);  
    
    % Compute error rotations in norm L2
    er_r = error_l2_norm_rot(rot, mc, element, coordinates, lambda);
    
    
    % Compute error div(sigma) in norm L2
    er_div = error_l2_norm_div(stress, mc, element, coordinates, lambda, G);
    
    % Print results
    fprintf(f, '%2.0f \t %6.5e \n', nelem, er_u);
    fprintf(ff, '%2.0f \t %6.5e \n', nelem, er_r);
    fprintf(fff, '%2.0f \t %6.5e \n', nelem, er_div);
    
end
fclose(f);
fclose(ff);
fclose(fff);

% Compute deformate
%def = defomesh(spost,element,coordinates) ;
% Plot solution
%plotsol(coordinates,element,stress,spost,rot,def,ndx,ndy) ;
% ------------------------------------------------------------------------ %
