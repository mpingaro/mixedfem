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
    
    %A = 2/(1+l);
    %B = 0.5*A*sin(pi*x)*sin(pi*y);
    %b = 1/25;

    %g(1,1) = b*(pi^2*(4*sin(2*pi*y)*(-1+2*cos(2*pi*x))-cos(pi*(x+y))+A*sin(pi*x)*sin(pi*y))); 
    %g(1,2) = b*(pi^2*(4*sin(2*pi*x)*( 1-2*cos(2*pi*y))-cos(pi*(x+y))+A*sin(pi*x)*sin(pi*y))); 

    g(1,1) = -pi^2*cos(pi*x)*sin(pi*y)*(l + mu + 2*l*cos(pi*y) + 12*mu*cos(pi*y));
    g(1,2) = -pi^2*sin(pi*x)*(l*cos(pi*y) + 3*mu*cos(pi*y) + 2*l*(2*cos(pi*y)^2 - 1) + 2*mu*(2*cos(pi*y)^2 - 1));   
    
    %g(1,1) = -4*mu*(2*y - 1)*(3*x^4 - 6*x^3 + 6*x^2*y^2 - 6*x^2*y + 3*x^2 - 6*x*y^2 + 6*x*y + y^2 - y);
    %g(1,2) = 4*mu*(2*x - 1)*(6*x^2*y^2 - 6*x^2*y + x^2 - 6*x*y^2 + 6*x*y - x + 3*y^4 - 6*y^3 + 3*y^2); 
    

    %% LOAD
    b_load(1,1) = b_load(1,1) + w*g(1,1)*DJ ;
    b_load(2,1) = b_load(2,1) + w*g(1,2)*DJ ;
end

end