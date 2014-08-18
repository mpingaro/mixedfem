function [gauss_p,gauss_w,npg] = gauss_quad()  

%% Gauss Quadrature
% - 9 Point Quadrature

% Direction x
gauss_p(1,1) = -0.774596669241483 ; 
gauss_p(1,2) = -0.774596669241483 ;
gauss_p(1,3) = -0.774596669241483 ; 
gauss_p(1,4) = 0.0 ; 
gauss_p(1,5) = 0.0 ; 
gauss_p(1,6) = 0.0 ; 
gauss_p(1,7) = 0.774596669241483 ; 
gauss_p(1,8) = 0.774596669241483 ; 
gauss_p(1,9) = 0.774596669241483 ;
% Direction y
gauss_p(2,1) = -0.774596669241483 ; 
gauss_p(2,2) = 0.0 ; 
gauss_p(2,3) = 0.774596669241483 ; 
gauss_p(2,4) = -0.774596669241483 ; 
gauss_p(2,5) = 0 ; 
gauss_p(2,6) = 0.774596669241483 ;
gauss_p(2,7) = -0.774596669241483 ; 
gauss_p(2,8) = 0.0 ; 
gauss_p(2,9) = 0.774596669241483 ;
% Weight of Quadrature
gauss_w(1,1) = 0.308641975308642 ; 
gauss_w(1,2) = 0.493827160493827 ;
gauss_w(1,3) = 0.308641975308642 ; 
gauss_w(1,4) = 0.493827160493827 ;
gauss_w(1,5) = 0.790123456790123 ; 
gauss_w(1,6) = 0.493827160493827 ;
gauss_w(1,7) = 0.308641975308642 ; 
gauss_w(1,8) = 0.493827160493827 ;
gauss_w(1,9) = 0.308641975308642 ;

npg = 9 ;

end
