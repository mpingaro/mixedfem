% Created by Marco Pingaro

function [J,JJ,DJ] = jacobian_quad(point,x,y)

% Gradient of shape functions
grad(1,1:4) = [-(1-y), 1-y, 1+y, -(1+y)].*0.25 ; % deriv along first direction
grad(2,1:4) = [-(1-x), -(1+x), 1+x, 1-x].*0.25 ; % deriv along second direction
    
% Jacobian Matrix
J(1,1) = grad(1,1:4)*point(1:4,1) ; % x_u
J(1,2) = grad(1,1:4)*point(1:4,2) ; % y_u
J(2,1) = grad(2,1:4)*point(1:4,1) ; % x_v
J(2,2) = grad(2,1:4)*point(1:4,2) ; % y_v

% Determinant of Jacobian Matrix
DJ = J(1,1)*J(2,2)-J(1,2)*J(2,1) ;
% Inverse of Jacobian Matrix
JJ(1,1) = J(2,2)/DJ ;  
JJ(1,2) = -J(1,2)/DJ ;
JJ(2,1) = -J(2,1)/DJ ; 
JJ(2,2) = J(1,1)/DJ ;

return