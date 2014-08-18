function [gauss_p,gauss_w,npg] = gauss_quad()  

% Gauss 6 points
x(1) = -0.932469514203152; x(6) = -x(1);
x(2) = -0.661209386466265; x(5) = -x(2);
x(3) = -0.238619186083197; x(4) = -x(3);
w(1) = 0.171324492379170; w(6) = w(1);
w(2) = 0.360761573048139; w(5) = w(2);
w(3) = 0.467913934572691; w(4) = w(3);


npg = 36;
gauss_p = zeros(2,npg) ;
gauss_w = zeros(1,npg) ;

np = 1 ; 
for i =1:6
    for j =1:6
        gauss_p(1,np) = x(i) ; 
        gauss_p(2,np) = x(j) ;
        gauss_w(1,np) = w(i)*w(j) ;

        np = np +1 ;
    end
end

return

%% Gauss Quadrature
% - 9 Point Quadrature

% Direction x
%gauss_p(1,1) = -0.774596669241483 ; 
%gauss_p(1,2) = -0.774596669241483 ;
%gauss_p(1,3) = -0.774596669241483 ; 
%gauss_p(1,4) = 0.0 ; 
%gauss_p(1,5) = 0.0 ; 
%gauss_p(1,6) = 0.0 ; 
%gauss_p(1,7) = 0.774596669241483 ; 
%gauss_p(1,8) = 0.774596669241483 ; 
%gauss_p(1,9) = 0.774596669241483 ;
% Direction y
%gauss_p(2,1) = -0.774596669241483 ; 
%gauss_p(2,2) = 0.0 ; 
%gauss_p(2,3) = 0.774596669241483 ; 
%gauss_p(2,4) = -0.774596669241483 ; 
%gauss_p(2,5) = 0 ; 
%gauss_p(2,6) = 0.774596669241483 ;
%gauss_p(2,7) = -0.774596669241483 ; 
%gauss_p(2,8) = 0.0 ; 
%gauss_p(2,9) = 0.774596669241483 ;
% Weight of Quadrature
%gauss_w(1,1) = 0.308641975308642 ; 
%gauss_w(1,2) = 0.493827160493827 ;
%gauss_w(1,3) = 0.308641975308642 ; 
%gauss_w(1,4) = 0.493827160493827 ;
%gauss_w(1,5) = 0.790123456790123 ; 
%gauss_w(1,6) = 0.493827160493827 ;
%gauss_w(1,7) = 0.308641975308642 ; 
%gauss_w(1,8) = 0.493827160493827 ;
%gauss_w(1,9) = 0.308641975308642 ;

%npg = 9 ;

%end
