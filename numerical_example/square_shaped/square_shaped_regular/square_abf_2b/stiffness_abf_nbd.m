% by Marco Pingaro

function [AELEM,BELEM,CELEM] = stiffness_abf_nbd(point,s,lambda,G)  

% Quadrature
[gauss_w, gauss_p] = GaussQuad2D(4,4);
npg = size(gauss_w,1);

%% - Legame
C = [(2*G + lambda)/(4*G*(G + lambda)), 0, 0, -lambda/(4*G*(G + lambda));...
    0, 1/(2*G), 0, 0; 0, 0, 1/(2*G), 0;...
    -lambda/(4*G*(G + lambda)), 0,0, (2*G + lambda)/(4*G*(G + lambda))];

%% ELEMENTARY MATRIX A & B
AELEM = zeros(16,16) ; 
BELEM = zeros(16,6) ;
CELEM = zeros(16,4) ;

for k = 1:npg
    x = gauss_p(k,1) ; y = gauss_p(k,2) ; w = gauss_w(k,1) ;

    [J,JJ,DJ] = jacobian_quad(point,x,y);
      
    %% --- Stress
    sig(:,1) = J*[0; -0.5+0.5*y]/DJ ;                % Shape 1 RT0
    sig(:,2) = J*[ 0.5+0.5*x; 0]/DJ ;                % Shape 2 RT0
    sig(:,3) = J*[ 0; 0.5+0.5*y]/DJ ;                % Shape 3 RT0
    sig(:,4) = J*[-0.5+0.5*x; 0]/DJ ;                % Shape 4 RT0
    % Bouble function 1
    sig(:,5) = [2*x*(y^2-1), 2*y*(x^2-1)]*JJ ;
    % Bouble function 2 
    sig(:,6) = [(1-3*x^2-2*x*y)*(1-y^2), (1-3*y^2-2*x*y)*(1-x^2)]*JJ ;

    % Bolle aggiuntive ABF
    sig(:,7) = [1-x^2; 0] ; % No mapping?
    sig(:,8) = [0; 1-y^2] ; % No mapping?

    % Tensor functions 
    sigt = zeros(2,2,16) ;
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
    sigt(:,:,13) = [sig(1,7), sig(2,7); 0, 0] ;
    sigt(:,:,14) = [sig(1,8), sig(2,8); 0, 0] ;
    sigt(:,:,15) = [0, 0; sig(1,7), sig(2,7)] ;
    sigt(:,:,16) = [0, 0; sig(1,8), sig(2,8)] ;
    
    sigma = reshape(sigt,4,16) ;

    %% --- Divergence
    div(:,1)  = s*[0.5, 0] ; 
    div(:,2)  = s*[0, 0.5] ; 
    div(:,3)  = s*[0.5, 0] ; 
    div(:,4)  = s*[0, 0.5] ; 
    div(:,5)  = s*[0.5, 0] ; 
    div(:,6)  = s*[0, 0.5] ; 
    div(:,7)  = s*[0.5, 0] ; 
    div(:,8)  = s*[0, 0.5] ; 
    div(:,9)  = [0, 0] ; 
    div(:,10) = [0, 0] ; 
    div(:,11) = [0, 0] ; 
    div(:,12) = [0, 0] ; 
    div(:,13) = [-2*x, 0] ; 
    div(:,14) = [-2*y, 0] ; 
    div(:,15) = [0, -2*x] ; 
    div(:,16) = [0, -2*y] ; 

    %% --- Displacement
    sp(:,1) = [1, 0] ;
    sp(:,2) = [x, 0] ;
    sp(:,3) = [y, 0] ;
    sp(:,4) = [0, 1] ;
    sp(:,5) = [0, x] ;
    sp(:,6) = [0, y] ; 

    %% --- Rotation
    rot(1,1) = 0.25*(1-x)*(1-y) ;
    rot(1,2) = 0.25*(1+x)*(1-y) ;
    rot(1,3) = 0.25*(1+x)*(1+y) ;
    rot(1,4) = 0.25*(1-x)*(1+y) ;
    asym = [rot; -rot] ;
    
    %% AELEM
    AELEM = AELEM + (C*sigma)'*sigma.*w.*DJ ;
    %% BELEM
    BELEM = BELEM + div'*sp.*w ;
    %% CELEM
    CELEM = CELEM + sigma([3,2],:)'*asym .*w.*DJ ;

end

end
