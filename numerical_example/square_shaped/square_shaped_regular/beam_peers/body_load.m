% by Marco Pingaro

% The body load is fix for the test

function b_load = body_load(p,l,m)  

% Different quadrature
%[gauss_p, gauss_w, npg] = quadrature_9() ;
%[gauss_p, gauss_w, npg] = quadrature_16() ;
%[gauss_p, gauss_w, npg] = quadrature_25() ;

[gauss_w, gauss_p] = GaussQuad2D(9,9);
npg = size(gauss_w,1);

b_load = zeros(2,1) ;
for k = 1:npg
    %x = gauss_p(1,k) ; y = gauss_p(2,k) ; w = gauss_w(1,k) ;
    x = gauss_p(k,1) ; y = gauss_p(k,2) ; w = gauss_w(k,1) ;

    % Gradient of shape functions
    grad(1,1:4) = [-(1-y), 1-y, 1+y, -(1+y)].*0.25 ; % deriv along first direction
    grad(2,1:4) = [-(1-x), -(1+x), 1+x, 1-x].*0.25 ; % deriv along second direction
  
    % Jacobian Matrix
    J(1,1) = sum( grad(1,1:4)*p(1:4,1) ) ; % x_u
    J(1,2) = sum( grad(2,1:4)*p(1:4,1) ) ; % x_v
    J(2,1) = sum( grad(1,1:4)*p(1:4,2) ) ; % y_u
    J(2,2) = sum( grad(2,1:4)*p(1:4,2) ) ; % y_v
       
    % Determinant of Jacobian Matrix
    DJ = J(1,1)*J(2,2)-J(1,2)*J(2,1) ;

    A = 2/(1+l);
    %B = 0.5*A*sin(pi*x)*sin(pi*y);
    b = 1/25;

    g(1,1) = b*(pi^2*(4*sin(2*pi*y)*(-1+2*cos(2*pi*x))-cos(pi*(x+y))+A*sin(pi*x)*sin(pi*y))); 
    g(1,2) = b*(pi^2*(4*sin(2*pi*x)*( 1-2*cos(2*pi*y))-cos(pi*(x+y))+A*sin(pi*x)*sin(pi*y))); 

    %% LOAD
    b_load(1,1) = b_load(1,1) + w*g(1,1)*DJ ;
    b_load(2,1) = b_load(2,1) + w*g(1,2)*DJ ;
end

end
