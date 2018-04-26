% by Marco Pingaro

function [AELEM,BELEM,CELEM] = stiffness_peers_nbd_old(point,s,cf)  

% Different quadrature
%[gauss_p, gauss_w, npg] = quadrature_9() ;
%[gauss_p, gauss_w, npg] = quadrature_16() ;
%[gauss_p, gauss_w, npg] = quadrature_25() ;
[gauss_w, gauss_p] = GaussQuad2D(4,4);
npg = size(gauss_w,1);

%% ELEMENTARY MATRIX A & B
AELEM = zeros(12,12) ; 
BELEM = zeros(12,2) ;
CELEM = zeros(12,4) ;

for k = 1:npg
    %x = gauss_p(1,k) ; y = gauss_p(2,k) ; w = gauss_w(1,k) ;
    x = gauss_p(k,1) ; y = gauss_p(k,2) ; w = gauss_w(k,1) ;

    % Gradient of shape functions
    grad(1,1:4) = [-(1-y), 1-y, 1+y, -(1+y)].*0.25 ; % deriv along first direction
    grad(2,1:4) = [-(1-x), -(1+x), 1+x, 1-x].*0.25 ; % deriv along second direction
    
    % Jacobian Matrix
    J(1,1) = sum( grad(1,1:4)*point(1:4,1) ) ; % x_u
    J(1,2) = sum( grad(2,1:4)*point(1:4,1) ) ; % x_v
    J(2,1) = sum( grad(1,1:4)*point(1:4,2) ) ; % y_u
    J(2,2) = sum( grad(2,1:4)*point(1:4,2) ) ; % y_v
       
    % Determinant of Jacobian Matrix
    DJ = J(1,1)*J(2,2)-J(1,2)*J(2,1) ;
    % Inverse transpose of Jacobian Matrix
    JJ(1,1) = J(2,2)/DJ ;  
    JJ(1,2) = -J(2,1)/DJ ;
    JJ(2,1) = -J(1,2)/DJ ; 
    JJ(2,2) = J(1,1)/DJ ;
      
    %% --- Stress
    sig(:,1) = J*[0; -0.5+0.5*y]/DJ ;                % Shape 1 RT0
    sig(:,2) = J*[ 0.5+0.5*x; 0]/DJ ;                % Shape 2 RT0
    sig(:,3) = J*[ 0; 0.5+0.5*y]/DJ ;                % Shape 3 RT0
    sig(:,4) = J*[-0.5+0.5*x; 0]/DJ ;                % Shape 4 RT0

    % Bouble function 1 (gradient)
    sig(:,5) = JJ*[2*x*(y^2-1); 2*y*(x^2-1)] ;
    % Bouble function 2 (gradient)
    sig(:,6) = JJ*[(1-3*x^2-2*x*y)*(1-y^2); (1-3*y^2-2*x*y)*(1-x^2)] ;

    % Tensor functions  
    sigt = zeros(2,2,12) ;
    sigt(:,:,1)  = [sig(1,1), sig(2,1); 0, 0].*s ;
    sigt(:,:,2)  = [0, 0; sig(1,1), sig(2,1)].*s ;
    sigt(:,:,3)  = [sig(1,2), sig(2,2); 0, 0].*s ;
    sigt(:,:,4)  = [0, 0; sig(1,2), sig(2,2)].*s ;
    sigt(:,:,5)  = [sig(1,3), sig(2,3); 0, 0].*s ;
    sigt(:,:,6)  = [0, 0; sig(1,3), sig(2,3)].*s ;
    sigt(:,:,7)  = [sig(1,4), sig(2,4); 0, 0].*s ;
    sigt(:,:,8)  = [0, 0; sig(1,4), sig(2,4)].*s ;
    sigt(:,:,9)  = [sig(2,5), -sig(1,5); 0, 0] ;
    sigt(:,:,10) = [0, 0; sig(2,5), -sig(1,5)] ;
    sigt(:,:,11) = [sig(2,6), -sig(1,6); 0, 0] ;
    sigt(:,:,12) = [0, 0; sig(2,6), -sig(1,6)] ;

    %% --- Rotation
    rot(1,1) = 0.25*(1-x)*(1-y) ;
    rot(1,2) = 0.25*(1+x)*(1-y) ;
    rot(1,3) = 0.25*(1+x)*(1+y) ;
    rot(1,4) = 0.25*(1-x)*(1+y) ;

    %% AELEM 
    for i = 1:12
        for j = 1:12
            AELEM(i,j) = AELEM(i,j) + ...
            ( cf(1,1)*trace( sigt(:,:,i)'*sigt(:,:,j) ) +...
            cf(1,2)*trace(sigt(:,:,i))*trace(sigt(:,:,j)) )*w*DJ ; 
        end
    end
    %% CELEM
    for i = 1:12
        for j = 1:4
            CELEM(i,j) = CELEM(i,j) +...
            rot(1,j)*( sigt(1,2,i) - sigt(2,1,i) )*w*DJ ;
        end
    end

end
%% BELEM bilinear form b(u, div(tau)) is costant for all elements!!
BELEM([1 3 5 7],1) = s*[2, 2, 2, 2] ;
BELEM([2 4 6 8],2) = s*[2, 2, 2, 2] ;

end
