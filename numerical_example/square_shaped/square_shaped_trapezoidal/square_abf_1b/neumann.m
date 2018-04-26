% by Marco Pingaro

function bn = neumann(ndx,ndy,g,side)

nedge = ndx*(ndy+1) + ndy*(ndx+1) ; 

if side == 1
    bc = 1:ndx*2 ;
    bc = sort(bc) ;
    bn = zeros(2,size(bc,2)) ;
    
    for i = 1 : size(bc,2)/2
        bn(2,(i-1)*2+[1 2]) = (-1)^(i)*[g(1,2) g(1,1)];
    end
    bn(1,:) = bc;

elseif side == 2
    bc = 2*ndx*ndy+2*(ndx+1)*ndy+1:2*nedge ;
    bc = sort(bc) ;
    bn = zeros(2,size(bc,2));
    
    for i = 1 : size(bc,2)/2
        bn(2,(i-1)*2+[1 2]) = (-1)^(i+ndy)*[g(1,2) g(1,1)];
    end
    bn(1,:) = bc;
    
elseif side == 3
    bc = 2*ndx+1:4*ndx+2:2*ndx*ndy+2*(ndx+1)*(ndy-1)+1 ;
    bc = [bc, bc+1] ;
    bc = sort(bc) ;
    bn = zeros(2,size(bc,2));
    
    for i = 1 : size(bc,2)/2
        bn(2,(i-1)*2+[1 2]) = (-1)^(i)*[g(1,1) g(1,2)];
    end
    bn(1,:) = bc;
    
elseif side == 4
    bc = 4*ndx+1:4*ndx+2:2*(nedge-ndx) ;
    bc = [bc, bc+1] ;
    bc = sort(bc) ;
    bn = zeros(2,size(bc,2));

    for i = 1 : size(bc,2)/2
        bn(2,(i-1)*2+[1 2]) = (-1)^(i+ndx)*[g(1,1) g(1,2)];
    end
    bn(1,:) = bc;
    
else
    fprintf('Warning: index out of bounds');
end

return

%bn1 = [bc1; repmat([g(1,2),g(1,1)],1, rep1)] ;
%bn2 = [bc2; repmat([g(2,2),g(2,1)],1, rep2)] ;
%bn3 = [bc3; repmat(g(3,:),1, rep3)] ;
%bn4 = [bc4; repmat(g(4,:),1, rep4)] ;