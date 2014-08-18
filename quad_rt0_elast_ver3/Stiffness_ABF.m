function [AELEM,BELEM,CELEM,b_load] = Stiffness_ABF(point,f,s,cf)  

%% Gauss Quadrature
% - 9 Point Quadrature

% Direction x
gauss_x(1,1) = -0.774596669241483 ; 
gauss_x(1,2) = -0.774596669241483 ;
gauss_x(1,3) = -0.774596669241483 ; 
gauss_x(1,4) = 0.0 ; 
gauss_x(1,5) = 0.0 ; 
gauss_x(1,6) = 0.0 ; 
gauss_x(1,7) = 0.774596669241483 ; 
gauss_x(1,8) = 0.774596669241483 ; 
gauss_x(1,9) = 0.774596669241483 ;
% Direction y
gauss_y(1,1) = -0.774596669241483 ; 
gauss_y(1,2) = 0.0 ; 
gauss_y(1,3) = 0.774596669241483 ; 
gauss_y(1,4) = -0.774596669241483 ; 
gauss_y(1,5) = 0 ; 
gauss_y(1,6) = 0.774596669241483 ;
gauss_y(1,7) = -0.774596669241483 ; 
gauss_y(1,8) = 0.0 ; 
gauss_y(1,9) = 0.774596669241483 ;
% Weight of Quadrature
w(1,1) = 0.308641975308642 ; 
w(1,2) = 0.493827160493827 ;
w(1,3) = 0.308641975308642 ; 
w(1,4) = 0.493827160493827 ;
w(1,5) = 0.790123456790123 ; 
w(1,6) = 0.493827160493827 ;
w(1,7) = 0.308641975308642 ; 
w(1,8) = 0.493827160493827 ;
w(1,9) = 0.308641975308642 ;

npg = 9 ;

%% ELEMENTARY MATRIX A & B
AELEM = zeros(10,10) ; b_load = zeros(2,1);
BELEM = zeros(10,2) ;
CELEM = zeros(10,4) ;
b_load = zeros(2,1) ;

for k = 1:npg
    x = gauss_x(1,k) ; 
    y = gauss_y(1,k) ;
    
    % Gradient of shape functions
    grad(1,1:4) = [-(1-y), 1-y, 1+y, -(1+y)].*0.25 ; % deriv along first direction
    grad(2,1:4) = [-(1-x), -(1+x), 1+x, 1-x].*0.25 ; % deriv along second direction
    
    % Jacobian Matrix
    J(1,1) = sum( grad(1,1:4)*point(1:4,1) ) ; % x_u
    J(1,2) = sum( grad(1,1:4)*point(1:4,2) ) ; % y_u
    J(2,1) = sum( grad(2,1:4)*point(1:4,1) ) ; % x_v
    J(2,2) = sum( grad(2,1:4)*point(1:4,2) ) ; % y_v
       
    % Determinant of Jacobian Matrix
    DJ = J(1,1)*J(2,2)-J(1,2)*J(2,1) ;
    % Inverse of Jacobian Matrix
    JJ(1,1) = J(2,2)/DJ ;  
    JJ(1,2) = -J(1,2)/DJ ;
    JJ(2,1) = -J(2,1)/DJ ; 
    JJ(2,2) = J(1,1)/DJ ;
      
    %% --- Stress
    sig(:,1) = J*[0; -0.5+0.5*y]/DJ ;                % Shape 1
    sig(:,2) = J*[0.5+0.5*x; 0]/DJ ;                 % Shape 2
    sig(:,3) = J*[0; 0.5+0.5*y]/DJ ;                 % Shape 3
    sig(:,4) = J*[-0.5+0.5*x; 0]/DJ ;                % Shape 4

    %sig(:,5) = J*[x*x-1; 0]/DJ;                     % Shape 5 Stabilization
    %sig(:,6) = J*[0; y*y-1]/DJ;                     % Shape 6 Stabilization

    sig(:,5) = JJ'*[2*x*(y^2-1); 2*y*(x^2-1)] ;      % Bouble function      
    sig(:,5) = [sig(2,5); -sig(1,5)]*DJ ;

    %% --- Rotation
    rot(1,1) =  0.25*(1-x)*(1-y) ;
    rot(1,2) =  0.25*(1+x)*(1-y) ;
    rot(1,3) =  0.25*(1+x)*(1+y) ;
    rot(1,4) =  0.25*(1-x)*(1+y) ;

    %% -- Displacement
    %disp(1) = [1; 0];
    %disp(2) = [0; 1];

    %% AELEM bilinear form a(C^(-1)*sigma,tau)
    % Upper triangular matrix 
    AELEM(1,1) = AELEM(1,1) + w(1,k)*(cf(1,1)*(sig(1,1)*sig(1,1)...
        +sig(2,1)*sig(2,1))+cf(1,2)*sig(1,1)*sig(1,1))*DJ ;
    
    AELEM(1,2) = AELEM(1,2) + w(1,k)*cf(1,2)*sig(1,1)*sig(2,1)*DJ ;
    
    AELEM(1,3) = AELEM(1,3) + w(1,k)*(cf(1,1)*(sig(1,1)*sig(1,2)...
        +sig(2,1)*sig(2,2))+cf(1,2)*sig(1,1)+sig(1,2))*DJ ;
    
    AELEM(1,4) = AELEM(1,4) + w(1,k)*cf(1,2)*sig(1,1)*sig(2,2)*DJ ;   
    
    AELEM(1,5) = AELEM(1,5) + w(1,k)*(cf(1,1)*(sig(1,1)*sig(1,3)...
        +sig(2,1)*sig(2,3))+cf(1,2)*sig(1,1)*sig(1,3))*DJ ;
    
    AELEM(1,6) = AELEM(1,6) + w(1,k)*cf(1,2)*sig(1,1)*sig(2,3)*DJ ;
    
    AELEM(1,7) = AELEM(1,7) + w(1,k)*(cf(1,1)*(sig(1,1)*sig(1,4)...
        +sig(2,1)*sig(2,4))+cf(1,2)*sig(1,1)*sig(1,4))*DJ ;
    
    AELEM(1,8) = AELEM(1,8) + w(1,k)*cf(1,2)*sig(1,1)*sig(2,4)*DJ ;
    
    AELEM(1,9) = AELEM(1,9) + w(1,k)*s*(cf(1,1)*(sig(1,1)*sig(1,5)...
        +sig(2,1)*sig(2,5))+cf(1,2)*sig(1,1)*sig(1,5))*DJ ;
    
    AELEM(1,10) = AELEM(1,10) + w(1,k)*s*cf(1,2)*sig(1,1)*sig(2,5)*DJ ;
    
    %AELEM(1,11) = AELEM(1,11) + w(1,k)*s*(cf(1,1)*(sig(1,1)*sig(1,6)...
    %    +sig(2,1)*sig(2,6))+cf(1,2)*sig(1,1)*sig(1,6))*DJ;
    
    %AELEM(1,12) = AELEM(1,12) + w(1,k)*s*cf(1,2)*sig(1,1)*sig(2,6)*DJ;
    
    %
    AELEM(2,2) = AELEM(2,2) + w(1,k)*(cf(1,1)*(sig(1,1)*sig(1,1)...
        +sig(2,1)*sig(2,1))+cf(1,2)*sig(2,1)*sig(2,1))*DJ ;
    
    AELEM(2,3) = AELEM(2,3) + w(1,k)*cf(1,2)*sig(2,1)*sig(1,2)*DJ ;
    
    AELEM(2,4) = AELEM(2,4) + w(1,k)*(cf(1,1)*(sig(1,1)*sig(1,2)...
        +sig(2,1)*sig(2,2))+cf(1,2)*sig(2,1)*sig(2,2))*DJ ;
    
    AELEM(2,5) = AELEM(2,5) + w(1,k)*cf(1,2)*sig(2,1)*sig(1,3)*DJ ;
    
    AELEM(2,6) = AELEM(2,6) + w(1,k)*(cf(1,1)*(sig(1,1)*sig(1,3)...
        +sig(2,1)*sig(2,3))+cf(1,2)*sig(2,1)*sig(2,3))*DJ ; 
    
    AELEM(2,7) = AELEM(2,7) + w(1,k)*cf(1,2)*sig(2,1)*sig(1,4)*DJ ;
    
    AELEM(2,8) = AELEM(2,8) + w(1,k)*(cf(1,1)*(sig(1,1)*sig(1,4)...
        +sig(2,1)*sig(2,4))+cf(1,2)*sig(2,1)*sig(2,4))*DJ ;
    
    AELEM(2,9) = AELEM(2,9) + w(1,k)*s*cf(1,2)*sig(2,1)*sig(1,5)*DJ ;
    
    AELEM(2,10) = AELEM(2,10) + w(1,k)*s*(cf(1,1)*(sig(1,1)*sig(1,5)...
        +sig(2,1)*sig(2,5))+cf(1,2)*sig(2,1)*sig(2,5))*DJ ;
    
    %AELEM(2,11) = AELEM(2,11) + w(1,k)*s*cf(1,2)*sig(2,1)*sig(1,6)*DJ;
    
    %AELEM(2,12) = AELEM(2,12) + w(1,k)*s*(cf(1,1)*(sig(1,1)*sig(1,6)...
    %    +sig(2,1)*sig(2,6))+cf(1,2)*sig(2,1)*sig(2,6))*DJ;
    %
    AELEM(3,3) = AELEM(3,3) + w(1,k)*(cf(1,1)*(sig(1,2)*sig(1,2)...
        +sig(2,2)*sig(2,2))+cf(1,2)*sig(1,2)*sig(1,2))*DJ;
    
    AELEM(3,4) = AELEM(3,4) + w(1,k)*cf(1,2)*sig(1,2)*sig(2,2)*DJ;
    
    AELEM(3,5) = AELEM(3,5) + w(1,k)*(cf(1,1)*(sig(1,2)*sig(1,3)...
        +sig(2,2)*sig(2,3))+cf(1,2)*sig(1,2)*sig(1,3))*DJ;
    
    AELEM(3,6) = AELEM(3,6) + w(1,k)*cf(1,2)*sig(1,2)*sig(2,3)*DJ;
    
    AELEM(3,7) = AELEM(3,7) + w(1,k)*(cf(1,1)*(sig(1,2)*sig(1,4)...
        +sig(2,2)*sig(2,4))+cf(1,2)*sig(1,2)*sig(1,4))*DJ;
    
    AELEM(3,8) = AELEM(3,8) + w(1,k)*cf(1,2)*sig(1,2)*sig(2,4)*DJ;
    
    AELEM(3,9) = AELEM(3,9) + w(1,k)*s*(cf(1,1)*(sig(1,2)*sig(1,5)...
        +sig(2,2)*sig(2,5))+cf(1,2)*sig(1,2)*sig(1,5))*DJ;
    
    AELEM(3,10) = AELEM(3,10) + w(1,k)*s*cf(1,2)*sig(1,2)*sig(2,5)*DJ;
    
    %AELEM(3,11) = AELEM(3,11) + w(1,k)*s*(cf(1,1)*(sig(1,2)*sig(1,6)...
    %    +sig(2,2)*sig(2,6))+cf(1,2)*sig(1,2)*sig(1,6))*DJ;
    
    %AELEM(3,12) = AELEM(3,12) + w(1,k)*s*cf(1,2)*sig(1,2)*sig(2,6)*DJ;
    %
    AELEM(4,4) = AELEM(4,4) + w(1,k)*(cf(1,1)*(sig(1,2)*sig(1,2)...
        +sig(2,2)*sig(2,2))+cf(1,2)*sig(2,2)*sig(2,2))*DJ;
    
    AELEM(4,5) = AELEM(4,5) + w(1,k)*cf(1,2)*sig(2,2)*sig(1,3)*DJ;
    
    AELEM(4,6) = AELEM(4,6) + w(1,k)*(cf(1,1)*(sig(1,2)*sig(1,3)...
        +sig(2,2)*sig(2,3))+cf(1,2)*sig(2,2)*sig(2,3))*DJ;
    
    AELEM(4,7) = AELEM(4,7) + w(1,k)*cf(1,2)*sig(2,2)*sig(1,4)*DJ;
    
    AELEM(4,8) = AELEM(4,8) + w(1,k)*(cf(1,1)*(sig(1,2)*sig(1,4)...
        +sig(2,2)*sig(2,4))+cf(1,2)*sig(2,2)*sig(2,4))*DJ;
    
    AELEM(4,9) = AELEM(4,9) + w(1,k)*s*cf(1,2)*sig(2,2)*sig(1,5)*DJ;
    
    AELEM(4,10) = AELEM(4,10) + w(1,k)*s*(cf(1,1)*(sig(1,2)*sig(1,5)...
        +sig(2,2)*sig(2,5))+cf(1,2)*sig(2,2)*sig(2,5))*DJ;
    
    %AELEM(4,11) = AELEM(4,11) + w(1,k)*s*cf(1,2)*sig(2,2)*sig(1,6)*DJ;
    
    %AELEM(4,12) = AELEM(4,12) + w(1,k)*s*(cf(1,1)*(sig(1,2)*sig(1,6)...
    %    +sig(2,2)*sig(2,6))+cf(1,2)*sig(2,2)*sig(2,6))*DJ;
    %
    AELEM(5,5) = AELEM(5,5) + w(1,k)*(cf(1,1)*(sig(1,3)*sig(1,3)...
        +sig(2,3)*sig(2,3))+cf(1,2)*sig(1,3)*sig(1,3))*DJ;
    
    AELEM(5,6) = AELEM(5,6) + w(1,k)*cf(1,2)*sig(1,3)*sig(2,3)*DJ;
    
    AELEM(5,7) = AELEM(5,7) + w(1,k)*(cf(1,1)*(sig(1,3)*sig(1,4)...
        +sig(2,3)*sig(2,4))+cf(1,2)*sig(1,3)*sig(1,4))*DJ;
    
    AELEM(5,8) = AELEM(5,8) + w(1,k)*cf(1,2)*sig(1,3)*sig(2,4)*DJ;
    
    AELEM(5,9) = AELEM(5,9) + w(1,k)*s*(cf(1,1)*(sig(1,3)*sig(1,5)...
        +sig(2,3)*sig(2,5))+cf(1,2)*sig(1,3)*sig(1,5))*DJ;
    
    AELEM(5,10) = AELEM(5,10) + w(1,k)*s*cf(1,2)*sig(1,3)*sig(2,5)*DJ;
    
    %AELEM(5,11) = AELEM(5,11) + w(1,k)*s*(cf(1,1)*(sig(1,3)*sig(1,6)...
    %    +sig(2,3)*sig(2,6))+cf(1,2)*sig(1,3)*sig(1,6))*DJ;
    
    %AELEM(5,12) = AELEM(5,12) + w(1,k)*s*cf(1,2)*sig(1,3)*sig(2,6)*DJ;
    %
    AELEM(6,6) = AELEM(6,6) + w(1,k)*(cf(1,1)*(sig(1,3)*sig(1,3)...
        +sig(2,3)*sig(2,3))+cf(1,2)*sig(2,3)*sig(2,3))*DJ;
    
    AELEM(6,7) = AELEM(6,7) + w(1,k)*cf(1,2)*sig(2,3)*sig(1,4)*DJ;
    
    AELEM(6,8) = AELEM(6,8) + w(1,k)*(cf(1,1)*(sig(1,3)*sig(1,4)...
        +sig(2,3)*sig(2,4))+cf(1,2)*sig(2,3)*sig(2,4))*DJ;
    
    AELEM(6,9) = AELEM(6,9) + w(1,k)*s*cf(1,2)*sig(2,3)*sig(1,5)*DJ;
    
    AELEM(6,10) = AELEM(6,10) + w(1,k)*s*(cf(1,1)*(sig(1,3)*sig(1,5)...
        +sig(2,3)*sig(2,5))+cf(1,2)*sig(2,3)*sig(2,5))*DJ;
    
    %AELEM(6,11) = AELEM(6,11) + w(1,k)*s*cf(1,2)*sig(2,3)*sig(1,6)*DJ;
    
    %AELEM(6,12) = AELEM(6,12) + w(1,k)*s*(cf(1,1)*(sig(1,3)*sig(1,6)...
    %    +sig(2,3)*sig(2,6))+cf(1,2)*sig(2,3)*sig(2,6))*DJ;
    %
    AELEM(7,7) = AELEM(7,7) + w(1,k)*(cf(1,1)*(sig(1,4)*sig(1,4)...
        +sig(2,4)*sig(2,4))+cf(1,2)*sig(1,4)*sig(1,4))*DJ;
    
    AELEM(7,8) = AELEM(7,8) + w(1,k)*cf(1,2)*sig(1,4)*sig(2,4)*DJ;
    
    AELEM(7,9) = AELEM(7,9) + w(1,k)*s*(cf(1,1)*(sig(1,4)*sig(1,5)...
        +sig(2,4)*sig(2,5))+cf(1,2)*sig(1,4)*sig(1,5))*DJ;
    
    AELEM(7,10) = AELEM(7,10) + w(1,k)*s*cf(1,2)*sig(1,4)*sig(2,5)*DJ;
    
    %AELEM(7,11) = AELEM(7,11) + w(1,k)*s*(cf(1,1)*(sig(1,4)*sig(1,6)...
    %    +sig(2,4)*sig(2,6))+cf(1,2)*sig(1,4)*sig(1,6))*DJ;
    
    %AELEM(7,12) = AELEM(7,12) + w(1,k)*s*cf(1,2)*sig(1,4)*sig(2,6)*DJ;
    %
    AELEM(8,8) = AELEM(8,8) + w(1,k)*(cf(1,1)*(sig(1,4)*sig(1,4)...
        +sig(2,4)*sig(2,4))+cf(1,2)*sig(2,4)*sig(2,4))*DJ;
    
    AELEM(8,9) = AELEM(8,9) + w(1,k)*s*cf(1,2)*sig(2,4)*sig(1,5)*DJ;
    
    AELEM(8,10) = AELEM(8,10) + w(1,k)*s*(cf(1,1)*(sig(1,4)*sig(1,5)...
        +sig(2,4)*sig(2,5))+cf(1,2)*sig(2,4)*sig(2,5))*DJ;
    
    %AELEM(8,11) = AELEM(8,11) + w(1,k)*s*cf(1,2)*sig(2,4)*sig(1,6)*DJ;
    
    %AELEM(8,12) = AELEM(8,12) + w(1,k)*s*(cf(1,1)*(sig(1,4)*sig(1,6)...
    %    +sig(2,4)*sig(2,6))+cf(1,2)*sig(2,4)*sig(2,6))*DJ;
    %
    AELEM(9,9) = AELEM(9,9) + w(1,k)*(cf(1,1)*(sig(1,5)*sig(1,5)...
        +sig(2,5)*sig(2,5))+cf(1,2)*sig(1,5)*sig(1,5))*DJ;
    
    AELEM(9,10) = AELEM(9,10) + w(1,k)*cf(1,2)*sig(1,5)*sig(2,5)*DJ;
    
    %AELEM(9,11) = AELEM(9,11) + w(1,k)*(cf(1,1)*(sig(1,5)*sig(1,6)...
    %    +sig(2,5)*sig(2,6))+cf(1,2)*sig(1,5)*sig(1,6))*DJ;
    
    %AELEM(9,12) = AELEM(9,12) + w(1,k)*cf(1,2)*sig(1,5)*sig(2,6)*DJ;
    %
    AELEM(10,10) = AELEM(10,10) + w(1,k)*(cf(1,1)*(sig(1,5)*sig(1,5)...
        +sig(2,5)*sig(2,5))+cf(1,2)*sig(2,5)*sig(2,5))*DJ;
    
    %AELEM(10,11) = AELEM(10,11) + w(1,k)*cf(1,2)*sig(2,5)*sig(1,6)*DJ;
    
    %AELEM(10,12) = AELEM(10,12) + w(1,k)*(cf(1,1)*(sig(1,5)*sig(1,6)...
    %    +sig(2,5)*sig(2,6))+cf(1,2)*sig(2,5)*sig(2,6))*DJ;
    %
    %AELEM(11,11) = AELEM(11,11) + w(1,k)*(cf(1,1)*(sig(1,6)*sig(1,6)...
    %    +sig(2,6)*sig(2,6))+cf(1,2)*sig(1,6)*sig(1,6))*DJ;
    
    %AELEM(11,12) = AELEM(11,12) + w(1,k)*cf(1,2)*sig(1,6)*sig(2,6)*DJ;
    % 
    %AELEM(12,12) = AELEM(12,12) + w(1,k)*(cf(1,1)*(sig(1,6)*sig(1,6)...
    %    +sig(2,6)*sig(2,6))+cf(1,2)*sig(2,6)*sig(2,6))*DJ;
    %
    %% CELEM bilinear form c(eta,as(tau))
    
    CELEM(1,1) = CELEM(1,1) + w(1,k)*s*(-rot(1,1)*sig(2,1))*DJ;
    CELEM(1,2) = CELEM(1,2) + w(1,k)*s*(-rot(1,2)*sig(2,1))*DJ;       
    CELEM(1,3) = CELEM(1,3) + w(1,k)*s*(-rot(1,3)*sig(2,1))*DJ;
    CELEM(1,4) = CELEM(1,4) + w(1,k)*s*(-rot(1,4)*sig(2,1))*DJ;
    %
    CELEM(2,1) = CELEM(2,1) + w(1,k)*s*( rot(1,1)*sig(1,1))*DJ;
    CELEM(2,2) = CELEM(2,2) + w(1,k)*s*( rot(1,2)*sig(1,1))*DJ;
    CELEM(2,3) = CELEM(2,3) + w(1,k)*s*( rot(1,3)*sig(1,1))*DJ;
    CELEM(2,4) = CELEM(2,4) + w(1,k)*s*( rot(1,4)*sig(1,1))*DJ;
    %
    CELEM(3,1) = CELEM(3,1) + w(1,k)*s*(-rot(1,1)*sig(2,2))*DJ; 
    CELEM(3,2) = CELEM(3,2) + w(1,k)*s*(-rot(1,2)*sig(2,2))*DJ;
    CELEM(3,3) = CELEM(3,3) + w(1,k)*s*(-rot(1,3)*sig(2,2))*DJ;
    CELEM(3,4) = CELEM(3,4) + w(1,k)*s*(-rot(1,4)*sig(2,2))*DJ;
    %
    CELEM(4,1) = CELEM(4,1) + w(1,k)*s*( rot(1,1)*sig(1,2))*DJ;
    CELEM(4,2) = CELEM(4,2) + w(1,k)*s*( rot(1,2)*sig(1,2))*DJ;
    CELEM(4,3) = CELEM(4,3) + w(1,k)*s*( rot(1,3)*sig(1,2))*DJ;
    CELEM(4,4) = CELEM(4,4) + w(1,k)*s*( rot(1,4)*sig(1,2))*DJ;
    %
    CELEM(5,1) = CELEM(5,1) + w(1,k)*s*(-rot(1,1)*sig(2,3))*DJ;
    CELEM(5,2) = CELEM(5,2) + w(1,k)*s*(-rot(1,2)*sig(2,3))*DJ;
    CELEM(5,3) = CELEM(5,3) + w(1,k)*s*(-rot(1,3)*sig(2,3))*DJ;
    CELEM(5,4) = CELEM(5,4) + w(1,k)*s*(-rot(1,4)*sig(2,3))*DJ;
    %
    CELEM(6,1) = CELEM(6,1) + w(1,k)*s*( rot(1,1)*sig(1,3))*DJ;
    CELEM(6,2) = CELEM(6,2) + w(1,k)*s*( rot(1,2)*sig(1,3))*DJ;
    CELEM(6,3) = CELEM(6,3) + w(1,k)*s*( rot(1,3)*sig(1,3))*DJ;
    CELEM(6,4) = CELEM(6,4) + w(1,k)*s*( rot(1,4)*sig(1,3))*DJ;
    %
    CELEM(7,1) = CELEM(7,1) + w(1,k)*s*(-rot(1,1)*sig(2,4))*DJ; 
    CELEM(7,2) = CELEM(7,2) + w(1,k)*s*(-rot(1,2)*sig(2,4))*DJ;
    CELEM(7,3) = CELEM(7,3) + w(1,k)*s*(-rot(1,3)*sig(2,4))*DJ;
    CELEM(7,4) = CELEM(7,4) + w(1,k)*s*(-rot(1,4)*sig(2,4))*DJ;
    %
    CELEM(8,1) = CELEM(8,1) + w(1,k)*s*( rot(1,1)*sig(1,4))*DJ;
    CELEM(8,2) = CELEM(8,2) + w(1,k)*s*( rot(1,2)*sig(1,4))*DJ;
    CELEM(8,3) = CELEM(8,3) + w(1,k)*s*( rot(1,3)*sig(1,4))*DJ;
    CELEM(8,4) = CELEM(8,4) + w(1,k)*s*( rot(1,4)*sig(1,4))*DJ;
    %
    CELEM(9,1) = CELEM(9,1) + w(1,k)*(-rot(1,1)*sig(2,5))*DJ;
    CELEM(9,2) = CELEM(9,2) + w(1,k)*(-rot(1,2)*sig(2,5))*DJ;
    CELEM(9,3) = CELEM(9,3) + w(1,k)*(-rot(1,3)*sig(2,5))*DJ;
    CELEM(9,4) = CELEM(9,4) + w(1,k)*(-rot(1,4)*sig(2,5))*DJ;
    %
    CELEM(10,1) = CELEM(10,1) + w(1,k)*( rot(1,1)*sig(1,5))*DJ;
    CELEM(10,2) = CELEM(10,2) + w(1,k)*( rot(1,2)*sig(1,5))*DJ;
    CELEM(10,3) = CELEM(10,3) + w(1,k)*( rot(1,3)*sig(1,5))*DJ;
    CELEM(10,4) = CELEM(10,4) + w(1,k)*( rot(1,4)*sig(1,5))*DJ;
    %
    %CELEM(11,1) = CELEM(11,1) + w(1,k)*(-rot(1,1)*sig(2,6))*DJ;
    %CELEM(11,2) = CELEM(11,2) + w(1,k)*(-rot(1,2)*sig(2,6))*DJ;
    %CELEM(11,3) = CELEM(11,3) + w(1,k)*(-rot(1,3)*sig(2,6))*DJ;
    %CELEM(11,4) = CELEM(11,4) + w(1,k)*(-rot(1,4)*sig(2,6))*DJ;
    %
    %CELEM(12,1) = CELEM(12,1) + w(1,k)*( rot(1,1)*sig(1,6))*DJ;
    %CELEM(12,2) = CELEM(12,2) + w(1,k)*( rot(1,2)*sig(1,6))*DJ;
    %CELEM(12,3) = CELEM(12,3) + w(1,k)*( rot(1,3)*sig(1,6))*DJ;
    %CELEM(12,4) = CELEM(12,4) + w(1,k)*( rot(1,4)*sig(1,6))*DJ;
    
    %% BELEM
    %BELEM(9,1)  = BELEM(9,1) + w(1,k)*2*x; 
    %BELEM(11,1) = BELEM(11,1) + w(1,k)*2*y;
    
    %% LOAD
    b_load(1,1) = b_load(1,1) + w(1,k)*f(1,1)*DJ;
    b_load(2,1) = b_load(2,1) + w(1,k)*f(2,1)*DJ;

end
%% BELEM bilinear form b(u, div(tau)) is costant for all elements!!
BELEM([1 3 5 7],1) = s*[2, 2, 2, 2];
BELEM([2 4 6 8],2) = s*[2, 2, 2, 2];
%BELEM(10,2) = BELEM(9,1);
%BELEM(12,2) = BELEM(11,1);
keyboard
end
