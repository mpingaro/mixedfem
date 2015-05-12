% by Marco Pingaro

function [AELEM,BELEM,CELEM,b_load] = stiffness_abf(point,f,s,lambda,G)  

%[gauss_p, gauss_w, npg] = quadrature_9() ;
[gauss_p, gauss_w, npg] = quadrature_16() ;

%% - Legame
C = [(2*G + lambda)/(4*G*(G + lambda)), 0, 0, -lambda/(4*G*(G + lambda));...
    0, 1/(2*G), 0, 0; 0, 0, 1/(2*G), 0;...
    -lambda/(4*G*(G + lambda)), 0,0, (2*G + lambda)/(4*G*(G + lambda))];

%% ELEMENTARY MATRIX A & B
AELEM = zeros(16,16) ; 
BELEM = zeros(16,6) ;
CELEM = zeros(16,4) ;
b_load = zeros(6,1) ;

for k = 1:npg
    x = gauss_p(1,k) ; y = gauss_p(2,k) ; w = gauss_w(1,k) ;
    
    [J,JJ,DJ] = jacobian_quad(point,x,y);
      
    %% --- Stress
    sig(:,1) = J*[0; -0.5+0.5*y]/DJ ;                % Shape 1 RT0
    sig(:,2) = J*[ 0.5+0.5*x; 0]/DJ ;                % Shape 2 RT0
    sig(:,3) = J*[ 0; 0.5+0.5*y]/DJ ;                % Shape 3 RT0
    sig(:,4) = J*[-0.5+0.5*x; 0]/DJ ;                % Shape 4 RT0
    % Bouble function 1
    sig(:,5) = JJ*[2*x*(y^2-1); 2*y*(x^2-1)] ;
    % Bouble function 2 
    sig(:,6) = JJ*[(1-3*x^2-2*x*y)*(1-y^2); (1-3*y^2-2*x*y)*(1-x^2)] ;

    % Bolle aggiuntive ABF
    sig(:,7) = JJ*[1-x^2; 0] ;
    sig(:,8) = JJ*[0; 1-y^2] ;

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

    %% LOAD
    b_load(1,1) = b_load(1,1) + w*f(1,1)*DJ ;
    b_load(2,1) = b_load(2,1) + w*x*f(1,1)*DJ ;
    b_load(3,1) = b_load(3,1) + w*y*f(1,1)*DJ ;
    
    b_load(4,1) = b_load(4,1) + w*f(2,1)*DJ ;
    b_load(5,1) = b_load(5,1) + w*x*f(2,1)*DJ ;
    b_load(6,1) = b_load(6,1) + w*y*f(2,1)*DJ ;

end

end

% OLD PART OF FUNCTION
% %% AELEM 
%     for i = 1:16
%         for j = 1:16
%             AELEM(i,j) = AELEM(i,j) + ...
%             ( cf(1,1)*(sigt(1,1,i)*sigt(1,1,j)+sigt(1,2,i)*sigt(1,2,j)+sigt(2,1,i)*sigt(2,1,j)+sigt(2,2,i)*sigt(2,2,j)) +...
%             cf(1,2)*( sigt(1,1,i)+sigt(2,2,i) )*( sigt(1,1,j)+sigt(2,2,j) ) )*w*DJ ; 
%         end
%     end
% 
%     %% BELEM
%     for i=1:16
%         for j=1:6
%             BELEM(i,j) = BELEM(i,j) + w*sp(:,j)'*div(:,i) ;
%         end
%     end
% 
%     %% CELEM
%     for i = 1:16
%         for j = 1:4
%             CELEM(i,j) = CELEM(i,j) +...
%             rot(1,j)*( sigt(1,2,i) - sigt(2,1,i) )*w*DJ ;
%         end
%     end