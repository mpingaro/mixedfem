% by Marco Pingaro
%

function [bc1,bc2,bc3,bc4] = dirichlet(ndx,ndy)

nedge = ndx*(ndy+1) + ndy*(ndx+1) ; 

bc1 = 1:ndx*2 ;
bc2 = 2*ndx*ndy+2*(ndx+1)*ndy+1:2*nedge ;

bc3 = 2*ndx+1:4*ndx+2:2*ndx*ndy+2*(ndx+1)*(ndy-1)+1 ;
bc3 = [bc3, bc3+1] ;

bc4 = 4*ndx+1:4*ndx+2:2*(nedge-ndx) ;
bc4 = [bc4, bc4+1] ;

bc1 = sort(bc1) ;
bc2 = sort(bc2) ;
bc3 = sort(bc3) ;
bc4 = sort(bc4) ;

return
