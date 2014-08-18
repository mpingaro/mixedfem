

%% --------------------------------------------------------------------- %%
% ----------------------- MIXED ELEMENT PEERS QUAD ---------------------- %
%                      FOR LINEAR ELASTICITY PROBLEM                      %
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
% Sigma, Tau are tensor   (2x2)                                           %
% u, v       are vector   (2x1)                                           %
% eta, psi   are costant  (1x1)                                           %
% ----------------------------------------------------------------------- %

%% INPUT DATI 
clear all; close all; clc;
length  = 1;                     % lunghezza trave
heigth  = 1;                     % altezza trave
young   = 1;                     % modulo di Young
poisson = 0.0;                   % modulo di Poisson
ndx     = 4;                     % numero suddivisioni in x
ndy     = 4;                     % numero suddivisioni in y
f(1,1) = 1;                      % load distribiuted direction x
f(2,1) = 1;                      % load distribiuted direction y
bc = [] ;                        % Dirichlet boundary condirion
% ----------------------------------------------------------------------- %
lambda = young*poisson/((1+poisson)*(1-2*poisson)) ;
mu = young/(2*(1+poisson)) ;
cf(1,1) = 1/(2*mu) ;
cf(1,2) = -lambda/(4*mu*(mu+lambda)) ;

% Geometry
[coordinates,element,mc] = BeamRectangular(length,heigth,ndx,ndy) ;
nelem = size(element,1) ; 
nnod  = size(coordinates,1) ;
ngdls = max(max(mc)) ;
ngdd = 2*nelem ;
ngdr = nnod ;
ngdlt = ngdls + ngdd + ngdr ;

% Assembly global system
[K,load] = assembly(coordinates,element,mc,cf,f) ; 

% Solve linear system
[stress,spost,rot] = solve(K,load,bc,ngdls,ngdd,ngdr) ;

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
    spost_x = spost(1:2:ngdd-1) ;
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
    spost_y = spost(2:2:ngdd) ; 
    surf(x,y,spost_y(i,1)*ones(4,4))
    hold on
end
axis equal
view(60,80)

% Rotation
X = zeros(1,ndx+1) ;
Y = zeros(1,ndy+1) ;
Z_rot = zeros(ndx+1,ndy+1) ;
figure
for i = 1:nelem
    X = [1:ndx+1] ;
    Y = [1:ndy+1] ;
    Z_rot = reshape(rot,ndy+1,ndx+1) ;
    hold on
end
surf(X,Y,Z_rot)
axis equal
view(60,80)
