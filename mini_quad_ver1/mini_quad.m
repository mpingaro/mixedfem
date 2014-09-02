
%% --------------------------------------------------------------------- %%
% ------------------------ MIXED ELEMENT MINI QUAD ---------------------- %
% ---------------   ( by Marco Pingaro & Daniele Boffi )   -------------- %
% Phd Student IUSS Pavia,                                                 %
% mail  : marco.pingaro@iusspavia.it                                      %
% ----------------------------------------------------------------------- %
% The problem is : laplace(u) = f                                         %
% Dirichlet condition is u_D = 0 on boundary domain                       %
%                                                                         %
% The Weak Formulation is:                                                %
% a(epsilon(u),epsilon(v)) + b(p,div(u) ) = ( f,v )                       %
%                            b(q,div(p) ) = 0                             %
% ----------------------------------------------------------------------- %

%% INPUT DATI 
clear all; close all; clc;
length  = 5;                     % lunghezza trave
heigth  = 1;                     % altezza trave
young   = 1;                     % modulo di Young
ndx     = 5;                     % numero suddivisioni in x
ndy     = 1;                     % numero suddivisioni in y
f(1,1) = 0;                      % load distribiuted direction x
f(2,1) = -1.0;                   % load distribiuted direction y
cf = 50 ;
% ----------------------------------------------------------------------- %

% Geometry
[coordinates,element,bc1,bc2,bc3,bc4] = beam(length,heigth,ndx,ndy) ;
nelem = size(element,1) ; 
nnod  = size(coordinates,1) ;
ngdld = 2*nnod+3*nelem ;
ngdlp = nnod ;
ngdlt = ngdld + ngdlp ;

% Assembly global system
[K,load] = assembly(coordinates,element,cf,f) ; 

% Dirichlet boundary conditions
bc = [bc3,bc4] ;

% Solve linear system
[spost,pres] = solve(K,load,bc,ngdld,ngdlp) ;

% Plot
plotsol(coordinates,spost,pres,ndx,ndy,nnod) ;
