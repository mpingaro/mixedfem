% by Marco Pingaro

function [K,F] = assembly(coor,el,mc,cf,f)

nelem = size(el,1) ;
nnod = size(coor,1) ;
ngdls = max(max(mc)) ;

%% LAPLACE RT_0 QUAD ELEMENT
ngdlt = ngdls + 6*nelem + nnod ;
K = sparse(ngdlt,ngdlt) ;
A = sparse(ngdls,ngdls) ;
B = sparse(ngdls,6*nelem) ;
C = sparse(ngdls,nnod) ;
F = sparse(ngdlt,1) ;
for k = 1:nelem
    %% Physic coordiantes points of the elements
    p(1,[1 2]) = coor( el(k,1),[1 2]) ;
    p(2,[1 2]) = coor( el(k,2),[1 2]) ;
    p(3,[1 2]) = coor( el(k,3),[1 2]) ;
    p(4,[1 2]) = coor( el(k,4),[1 2]) ;
    s = mc(k,16) ;
    % Compute the bilinear form a(sigma,tau), b(u, div(tau)) & the load (f,v) 
    [AELEM,BELEM,CELEM,load] = stiffness_abf(p,f,s,cf) ;
    for i = 1:15
        for j = 1:15
            A(mc(k,i),mc(k,j)) = A(mc(k,i),mc(k,j)) + AELEM(i,j) ;
        end
        for jj=1:4
            C(mc(k,i),el(k,jj)) = C(mc(k,i),el(k,jj)) + CELEM(i,jj) ;
        end
        B(mc(k,i),6*(k-1)+[1, 2, 3, 4, 5, 6]) = BELEM(i,[1, 2, 3, 4, 5, 6]) ;
        F(ngdls+6*(k-1)+[1, 2, 3, 4, 5, 6],1) = -load([1, 2, 3, 4, 5, 6],1) ;
    end
end
% GLOBAL SYSTEM
K = [A, B, C; B',sparse(6*nelem,6*nelem+nnod); C', sparse(nnod,6*nelem+nnod)] ;

return
