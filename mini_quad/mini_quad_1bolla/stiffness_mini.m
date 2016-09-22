%by Marco Pingaro

function [AELEM,BELEM,b_load] = stiffness_mini(point,f,cf)  

[gauss_p, gauss_w, npg] = quadrature_9() ;

%% ELEMENTARY MATRIX A & B
AELEM = zeros(11,11) ; 
BELEM = zeros(11,4) ;
b_load = zeros(11,1) ;

for k = 1:npg
    x = gauss_p(1,k) ; 
    y = gauss_p(2,k) ;
    w = gauss_w(1,k) ;

    % Shape functions
    
    psi(1,1) = 0.25*(1-x)*(1-y) ;
    psi(1,2) = 0.25*(1+x)*(1-y) ;
    psi(1,3) = 0.25*(1+x)*(1+y) ;
    psi(1,4) = 0.25*(1-x)*(1+y) ;
    % Bouble functions 1 versione
    psi(1,5) = 4*psi(1,1)*(1-x^2)*(1-y^2) ; 

    % Gradient of shape functions 
    
    % deriv along first direction
    grad(1,1:5) = [-(1-y), 1-y, 1+y, -(1+y), 4*(-2*x-1+3*x^2)*(1-y-y^2+y^3)].*0.25 ;
    % deriv along second direction
    grad(2,1:5) = [-(1-x), -(1+x), 1+x, 1-x, 4*(-2*y-1+3*y^2)*(1-x-x^2+x^3)].*0.25 ;
     
    J(1,1) = sum( grad(1,1:4)*point(1:4,1) ) ; % x_u
    J(1,2) = sum( grad(2,1:4)*point(1:4,1) ) ; % x_v
    J(2,1) = sum( grad(1,1:4)*point(1:4,2) ) ; % y_u
    J(2,2) = sum( grad(2,1:4)*point(1:4,2) ) ; % y_v
  
    % Determinant of Jacobian Matrix
    DJ = J(1,1)*J(2,2)-J(1,2)*J(2,1) ;
    % Inverse transpose of Jacobian Matrix
    JJ(1,1) =  J(2,2)/DJ ;  
    JJ(1,2) = -J(2,1)/DJ ;
    JJ(2,1) = -J(1,2)/DJ ; 
    JJ(2,2) =  J(1,1)/DJ ;
    % Trasformation gradient  
    grad_u = JJ*grad ;

    % Strain
    epsi(:,:,1) = [grad_u(1,1), grad_u(2,1)*0.5; grad_u(2,1)*0.5, 0] ;
    epsi(:,:,2) = [0, grad_u(1,1)*0.5; grad_u(1,1)*0.5, grad_u(2,1)] ;
    epsi(:,:,3) = [grad_u(1,2), grad_u(2,2)*0.5; grad_u(2,2)*0.5, 0] ;
    epsi(:,:,4) = [0, grad_u(1,2)*0.5; grad_u(1,2)*0.5, grad_u(2,2)] ;
    epsi(:,:,5) = [grad_u(1,3), grad_u(2,3)*0.5; grad_u(2,3)*0.5, 0] ;
    epsi(:,:,6) = [0, grad_u(1,3)*0.5; grad_u(1,3)*0.5, grad_u(2,3)] ;
    epsi(:,:,7) = [grad_u(1,4), grad_u(2,4)*0.5; grad_u(2,4)*0.5, 0] ;
    epsi(:,:,8) = [0, grad_u(1,4)*0.5; grad_u(1,4)*0.5, grad_u(2,4)] ;
    epsi(:,:,9) = [grad_u(1,5), grad_u(2,5)*0.5; grad_u(2,5)*0.5, 0] ;
    epsi(:,:,10) = [0, grad_u(1,5)*0.5; grad_u(1,5)*0.5, grad_u(2,5)] ;

    %% AELEM 
    for i = 1:10
        for j = 1:10
            AELEM(i,j) = AELEM(i,j) + ...
            cf*(epsi(1,1,i)*epsi(1,1,j) + epsi(2,2,i)*epsi(2,2,j)+...
            2*epsi(1,2,i)*epsi(1,2,j) )*w*DJ ;
        end
    end

    %% LOAD
    b_load(1,1) = b_load(1,1) + w * f(1,1) * psi(1,1) * DJ ;
    b_load(2,1) = b_load(2,1) + w * f(2,1) * psi(1,1) * DJ ;
    b_load(3,1) = b_load(3,1) + w * f(1,1) * psi(1,2) * DJ ;
    b_load(4,1) = b_load(4,1) + w * f(2,1) * psi(1,2) * DJ ;
    b_load(5,1) = b_load(5,1) + w * f(1,1) * psi(1,3) * DJ ;
    b_load(6,1) = b_load(6,1) + w * f(2,1) * psi(1,3) * DJ ;
    b_load(7,1) = b_load(7,1) + w * f(1,1) * psi(1,4) * DJ ;
    b_load(8,1) = b_load(8,1) + w * f(2,1) * psi(1,4) * DJ ;
    b_load(9,1) = b_load(9,1) + w * f(1,1) * psi(1,5) * DJ ;
    b_load(10,1) = b_load(10,1) + w * f(2,1) * psi(1,5) * DJ ;

    %% BELEM
    for nb = 1:4
        BELEM(1,nb) = BELEM(1,nb) + w *( psi(1,nb)*grad(1,1) ) ;
        BELEM(2,nb) = BELEM(2,nb) + w *( psi(1,nb)*grad(2,1) ) ;
        BELEM(3,nb) = BELEM(3,nb) + w *( psi(1,nb)*grad(1,2) ) ;
        BELEM(4,nb) = BELEM(4,nb) + w *( psi(1,nb)*grad(2,2) ) ;
        BELEM(5,nb) = BELEM(5,nb) + w *( psi(1,nb)*grad(1,3) ) ;
        BELEM(6,nb) = BELEM(6,nb) + w *( psi(1,nb)*grad(2,3) ) ;
        BELEM(7,nb) = BELEM(7,nb) + w *( psi(1,nb)*grad(1,4) ) ;
        BELEM(8,nb) = BELEM(8,nb) + w *( psi(1,nb)*grad(2,4) ) ;
        BELEM(9,nb) = BELEM(9,nb) + w *( psi(1,nb)*grad(1,5) ) ;
        BELEM(10,nb) = BELEM(10,nb) + w *( psi(1,nb)*grad(2,5) ) ;
    end
end

return
