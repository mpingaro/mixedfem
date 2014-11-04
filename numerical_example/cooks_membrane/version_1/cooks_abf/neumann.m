% by Marco Pingaro

function [bn1,bn2,bn3,bn4] = neumann(ndx,ndy,g)

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

rep1 = size(bc1,2)/2 ;
rep2 = size(bc2,2)/2 ;
rep3 = size(bc3,2)/2 ;
rep4 = size(bc4,2)/2 ;
bn1 = [bc1; repmat([g(1,2), g(1,1)],1, rep1)] ;
bn2 = [bc2; repmat([g(2,2), g(2,1)],1, rep2)] ;
bn3 = [bc3; repmat(g(3,:),1, rep3)] ;
bn4 = [bc4; repmat(g(4,:),1, rep4)] ;

return
