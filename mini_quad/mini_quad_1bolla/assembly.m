% by Marco Pingaro

function [K,F] = assembly(coor,el,cf,f)

nelem = size(el,1) ;
nnod = size(coor,1) ;
nd = 2*nnod+2*nelem ;
nt = nd + nnod ;

%% LAPLACE RT_0 QUAD ELEMENT
K = sparse(nt,nt) ;
A = sparse(nd,nd) ;
B = sparse(nd,nnod) ;
F = sparse(nt,1) ;
mc = zeros(1,10) ;
for k = 1:nelem
    %% Physic coordiantes points of the elements
    p(1,[1 2]) = coor( el(k,1),[1 2]) ;
    p(2,[1 2]) = coor( el(k,2),[1 2]) ;
    p(3,[1 2]) = coor( el(k,3),[1 2]) ;
    p(4,[1 2]) = coor( el(k,4),[1 2]) ;
    % Compute the bilinear form a(epsilon,epsilon), b(p, div(u)) & the load (f,v) 
    [AELEM,BELEM,load] = stiffness_mini(p,f,cf) ;

    % Corrispondence matrix    
    mc([1,2]) = [2*el(k,1)-1, 2*el(k,1)] ;
    mc([3,4]) = [2*el(k,2)-1, 2*el(k,2)] ;
    mc([5,6]) = [2*el(k,3)-1, 2*el(k,3)] ;
    mc([7,8]) = [2*el(k,4)-1, 2*el(k,4)] ;
    ind = 2*nnod+2*(k-1) + 1 ;
    mc([9,10]) = [ind, ind+1];

    for i = 1:10
        for j = 1:10
            A(mc(1,i),mc(1,j)) = A(mc(1,i),mc(1,j)) + AELEM(i,j) ;
        end
        
        F(mc(1,i),1) = F(mc(1,i),1) + load(i,1) ;
        
        for jj=1:4
            B(mc(1,i),el(k,jj)) = B(mc(1,i),el(k,jj)) - BELEM(i,jj) ;
        end
        
    end
end
% GLOBAL SYSTEM
%
K = [A, B; B',sparse(nnod,nnod)] ;

return
