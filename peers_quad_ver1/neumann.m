% by Marco Pingaro
%

function [bc1,bc2,bc3,bc4] = neumann(ndx,ndy)

nedge = ndx*(ndy+1) + ndy*(ndx+1) ; 

bc1 = 1:ndx*2 ;
bc2 = 2*ndx*ndy+2*(ndx+1)*ndy+1:2*nedge ;
bc3 = [];
bc4 = [];

return
