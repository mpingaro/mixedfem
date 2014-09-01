

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
% a(C^(-1)Sigma,Tau ) + b( u,div(Tau) ) + c( eta, as(tau) ) = 0           %
%                          b( v, div(Sigma) )               = ( f,v )     %
%                          c( psi, as(sigma) )              = 0           %
% ----------------------------------------------------------------------- %

%% INPUT DATI 
clear all; close all; clc;
length  = 5;                     % lunghezza trave
heigth  = 1;                     % altezza trave
young   = 1;                     % modulo di Young
ndx     = 2;                     % numero suddivisioni in x
ndy     = 1;                     % numero suddivisioni in y
f(1,1) = 0;                      % load distribiuted direction x
f(2,1) = -0.1;                      % load distribiuted direction y
cf = 1 ;
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
bc = bc3 ;

% Solve linear system
[spost,pres] = solve(K,load,bc,ngdld,ngdlp) ;

% PLOT SOLUTION DISPLACEMENT   
x = reshape(coordinates(:,1),ndx+1,ndy+1) ;
y = reshape(coordinates(:,2),ndx+1,ndy+1) ;
figure, mesh(x,y,zeros(ndx+1,ndy+1))
axis equal
view(0,90)

% Displacement X
figure
for i = 1 :nelem
    x = coordinates(element(i,1:4),1) ;
    y = coordinates(element(i,1:4),2) ;
    %mesh(x,y,spost_x(i,1)*ones(4,4))
    spost_x = spost(1:2:2*nnod-1) ;
    surf(x,y,spost_x(i,1)*ones(4,4))
    hold on
end
axis equal
view(60,80)

% Displacement Y
figure
for i = 1:nelem
    x = coordinates(element(i,1:4),1) ;
    y = coordinates(element(i,1:4),2) ;
    %mesh(x,y,spost_y(i,1)*ones(4,4))
    spost_y = spost(2:2:2*nnod) ; 
    surf(x,y,spost_y(i,1)*ones(4,4))
    hold on
end
axis equal
view(60,80)

% Deformation
x = reshape(coordinates(:,1),ndx+1,ndy+1) ;
y = reshape(coordinates(:,2),ndx+1,ndy+1) ;
spost_x = reshape(spost_x,ndx+1,ndy+1) ;
spost_y = reshape(spost_y,ndx+1,ndy+1) ;
def_x = x+spost_x ;
def_y = y+spost_y ;
figure, mesh(def_x,def_y,zeros(ndx+1,ndy+1))
axis equal
view(0,90)

% Pression
X = zeros(1,ndx+1) ;
Y = zeros(1,ndy+1) ;
Z_rot = zeros(ndx+1,ndy+1) ;
figure
for i = 1:nelem
    X = [1:ndx+1] ;
    Y = [1:ndy+1] ;
    pres = reshape(pres,ndy+1,ndx+1) ;
    hold on
end
surf(X,Y,pres)
axis equal
view(60,80)
