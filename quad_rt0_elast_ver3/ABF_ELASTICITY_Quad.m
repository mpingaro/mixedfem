%% --------------------------------------------------------------------- %%
%                                                                         %
% ----------------------- MIXED ELEMENT ABF_0 QUAD ---------------------- %
%                      FOR LINEAR ELASTICITY PROBLEM                      %
% -----------------   ( by Dott. Ing. Marco Pingaro )   ----------------- %
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
%                                                                         %
%                                                                         %
% Shape Function                                                          %
%                                                                         %
% Stress: RT_0 = [a1 + b1*x , c1 + d1*y; a2 + b2*x , c2 + d2*y]           %
%              + [x^2-1 , 0; 0, 0] + [0, 0; x^2-1, 0]                     %
%              + [0 , y^2-1; 0, 0] + [0, 0; y^2-1, 0]                     %
% Displacement : P_0 = [a, 0] + [0, b]                      Discontinuous %
% Rotation     : P_1 = a*x + b*y + c*x*y + d                              %
%                                                                         %
% Element Reference Domain (-1, 1) x (-1, 1)                              %
% ----------------------------------------------------------------------- %

%% INPUT DATI 
clear all; close all; clc;
length  = 1;                     % lunghezza trave
heigth  = 1;                     % altezza trave
young   = 1;                     % modulo di Young
poisson = 0.0;                   % modulo di Poisson
ndx     = 3;                     % numero suddivisioni in x
ndy     = 3;                     % numero suddivisioni in y
f_load(1,1) = 1;                 % load distribiuted direction x
f_load(2,1) = 1;                 % load distribiuted direction y

%%
lambda = young*poisson/((1+poisson)*(1-2*poisson)) ;
mu = young/(2*(1+poisson)) ;
cf(1,1) = 1/(2*mu) ;
cf(1,2) = -lambda/(4*mu*(mu+lambda)) ;

%%
[coordinates,element,mc] = BeamRectangular(length,heigth,ndx,ndy) ;
nelem = size(element,1) ;
nnod = size(coordinates,1) ;
ngdls = max(max(mc)) ;

%% LAPLACE RT_0 QUAD ELEMENT
ngdlt = ngdls + 2*nelem + nnod ;
K = sparse(ngdlt,ngdlt) ;
A = sparse(ngdls,ngdls) ;
B = sparse(ngdls,2*nelem) ;
C = sparse(ngdls,nnod) ;
load = sparse(ngdlt,1) ;
for k = 1:nelem
    %% Physic coordiantes points of the elements
    point(1,[1 2]) = coordinates( element(k,1),[1 2]) ;
    point(2,[1 2]) = coordinates( element(k,2),[1 2]) ;
    point(3,[1 2]) = coordinates( element(k,3),[1 2]) ;
    point(4,[1 2]) = coordinates( element(k,4),[1 2]) ;
    segno = mc(k,13) ;
    % Compute the bilinear form a(sigma,tau), b(u, div(tau)) & the load (f,v) 
    [AELEM,BELEM,CELEM,b_load] = Stiffness_ABF_NEW(point,f_load,segno,cf) ;
    for i = 1:12
        for j = 1:12
            %if i>j
            %    AEL = AELEM(j,i) ;
            %else
            %    AEL = AELEM(i,j) ;
            %end
            A(mc(k,i),mc(k,j)) = A(mc(k,i),mc(k,j)) + AELEM(i,j) ;
        end
        for jj=1:4
            C(mc(k,1),element(k,jj)) = C(mc(k,1),element(k,jj)) + CELEM(i,jj) ;
        end

        for ii = 1:2
            B(mc(k,i),2*(k-1)+ii) = BELEM(i,ii) ;
            load(ngdls+2*(k-1)+ii,1) = b_load(ii,1) ;
        end
    end
end
% GLOBAL SYSTEM
K = [A, B, C; B',sparse(2*nelem,2*nelem+nnod); C', sparse(nnod,2*nelem+nnod)] ;

%% ASSIGN BOUNDARY CONDITION (sigma * n = 0)
I1 = 1:2*ndx ;
I2 = 2*ndx*ndy+2*ndy*(ndx+1)+1:2*ndx*(ndy+1)+2*ndy*(ndx+1) ;
I = [] ;
for jj = 1: size(I,2)
    K(I(1,jj),:) = 0 ;
    K(I(1,jj),I(1,jj)) = 1 ;
end

% SOLUTION OF LINEAR SYSTEM
soluz = K\load ;
spost = soluz(ngdls+1:ngdls+2*nelem) ;
rot = soluz(ngdls+2*nelem+1:ngdlt) ;

spost_x = spost(1:2:2*nelem-1) ;  % Costant part of displacement
spost_y = spost(2:2:2*nelem) ;    % Costant part of displacement

% PLOT SOLUTION DISPLACEMENT   
x = reshape(coordinates(:,1),ndx+1,ndy+1) ;
y = reshape(coordinates(:,2),ndx+1,ndy+1) ;
figure, mesh(x,y,zeros(ndx+1,ndy+1))
axis equal
view(0,90)

% Displacement X
figure
for i = 1:size(element,1)
    x = coordinates(element(i,1:4),1) ;
    y = coordinates(element(i,1:4),2) ;
    %mesh(x,y,spost_x(i,1)*ones(4,4))
    surf(x,y,spost_x(i,1)*ones(4,4))
    hold on
end
axis equal
view(60,80)

% Displacement Y
figure
for i = 1:size(element,1)
    x = coordinates(element(i,1:4),1) ;
    y = coordinates(element(i,1:4),2) ;
    %mesh(x,y,spost_y(i,1)*ones(4,4))
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
for i = 1:size(element,1)
    X = [1:ndx+1] ;
    Y = [1:ndy+1] ;
    Z_rot = reshape(rot,ndy+1,ndx+1) ;
    hold on
end
surf(X,Y,Z_rot)
axis equal
view(60,80)
