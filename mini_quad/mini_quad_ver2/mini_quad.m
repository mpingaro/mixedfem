
%% --------------------------------------------------------------------- %%
% ------------------------ MIXED ELEMENT MINI QUAD ---------------------- %
% -----------------------   ( by Marco Pingaro )   ---------------------- %
% Phd Student IUSS Pavia,                                                 %
% mail  : marco.pingaro@iusspavia.it                                      %
% ----------------------------------------------------------------------- %
% The problem is : laplace(u) = f                                         %
% Dirichlet condition is u_D = 0 on boundary domain                       %
%                                                                         %
% The Weak Formulation is:                                                %
% a(epsilon(u),epsilon(v)) - b(p,div(u) ) = ( f,v )                       %
%                          - b(q,div(p) ) = 0                             %
% ----------------------------------------------------------------------- %

%% INPUT DATI 
clear all; close all; clc;
length  = 4;                     % lunghezza trave
heigth  = 2;                     % altezza trave
young   = 50;                    % modulo di Young
ndx     = 2;                     % numero suddivisioni in x
ndy     = 1;                     % numero suddivisioni in y
f(1,1) =  0.00;                  % load distribiuted direction x
f(2,1) = -0.50;                  % load distribiuted direction y
% ----------------------------------------------------------------------- %

% Geometry
[coordinates,element,bc1,bc2,bc3,bc4] = beam(length,heigth,ndx,ndy) ;
nelem = size(element,1) ; 
nnod  = size(coordinates,1) ;
ngdld = 2*nnod+4*nelem ;
ngdlp = nnod ;
ngdlt = ngdld + ngdlp ;
cf = young/1.5 ;
% Assembly global system
[K,load] = assembly(coordinates,element,cf,f) ; 

% Dirichlet boundary conditions
bc = [bc3] ;

% Solve linear system
[spost,pres] = solve(K,load,bc,ngdld,ngdlp) ;

% Plot
%plotsol(coordinates,spost,pres,ndx,ndy,nnod) ;
