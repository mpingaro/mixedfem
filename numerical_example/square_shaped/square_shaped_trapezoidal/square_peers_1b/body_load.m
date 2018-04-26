% by Marco Pingaro

% The body load is fix for the test

function b_load = body_load(p,l,mu)  

% Quadrature
[gauss_w, gauss_p] = GaussQuad2D(4,4);
npg = size(gauss_w,1);

b_load = zeros(2,1) ;
for k = 1:npg
    xi = gauss_p(k,1) ; eta = gauss_p(k,2) ; w = gauss_w(k,1) ;

    [J,JJ,DJ] = jacobian_quad(p,xi,eta);
    [x,y] = map_quad(p,gauss_p(k,:));
    
    g(1,1) = -pi^2*cos(pi*x)*sin(pi*y)*(l + mu + 2*l*cos(pi*y) + 12*mu*cos(pi*y));
    g(1,2) = -pi^2*sin(pi*x)*(l*cos(pi*y) + 3*mu*cos(pi*y) + 2*l*(2*cos(pi*y)^2 - 1) + 2*mu*(2*cos(pi*y)^2 - 1));    
    
    %% LOAD
    b_load(1,1) = b_load(1,1) + w*g(1,1)*DJ ;
    b_load(2,1) = b_load(2,1) + w*g(1,2)*DJ ;
end

end