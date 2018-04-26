% by Marco Pingaro

function [AELEM,BELEM,CELEM] = stiffness_peers_nbd(point,s,lambda,G)  

% Quadrature
[gauss_w, gauss_p] = GaussQuad2D(4,4);
npg = size(gauss_w,1);

%% - Legame
C = [(2*G + lambda)/(4*G*(G + lambda)), 0, 0, -lambda/(4*G*(G + lambda));...
    0, 1/(2*G), 0, 0; 0, 0, 1/(2*G), 0;...
    -lambda/(4*G*(G + lambda)), 0,0, (2*G + lambda)/(4*G*(G + lambda))];

%% ELEMENTARY MATRIX A & B
AELEM = zeros(12,12) ; 
BELEM = zeros(12,2) ;
CELEM = zeros(12,4) ;

for k = 1:npg
    x = gauss_p(k,1) ; y = gauss_p(k,2) ; w = gauss_w(k,1) ;

    [J,JJ,DJ] = jacobian_quad(point,x,y);
    
    %% --- Stress
    sig(:,1) = J*[0; -0.5+0.5*y]/DJ ;                % Shape 1 RT0
    sig(:,2) = J*[ 0.5+0.5*x; 0]/DJ ;                % Shape 2 RT0
    sig(:,3) = J*[ 0; 0.5+0.5*y]/DJ ;                % Shape 3 RT0
    sig(:,4) = J*[-0.5+0.5*x; 0]/DJ ;                % Shape 4 RT0

    % Bouble function 1 (gradient)
    sig(:,5) = [2*x*(y^2-1), 2*y*(x^2-1)]*JJ ;
    % Bouble function 2 (gradient)
    sig(:,6) = [(1-3*x^2-2*x*y)*(1-y^2), (1-3*y^2-2*x*y)*(1-x^2)]*JJ ;

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
    
    sigma = reshape(sigt,4,12);
           
    %% AELEM
    AELEM = AELEM + (C*sigma)'*sigma.*w.*DJ ;
    %% CELEM
    asym = [rot; -rot];
    CELEM = CELEM + sigma([3,2],:)'*asym .*w.*DJ ;

end
%% BELEM bilinear form b(u, div(tau)) is costant for all elements!!
BELEM([1 3 5 7],1) = s*[2, 2, 2, 2] ;
BELEM([2 4 6 8],2) = s*[2, 2, 2, 2] ;

end