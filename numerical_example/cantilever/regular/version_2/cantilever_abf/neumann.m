% by Marco Pingaro

function [bn1,bn2,bn3,bn4] = neumann(cr,ndx,ndy,g,h)

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

for i = 1 : rep1
    bn1(2,(i-1)*2+[1 2]) = (-1)^(i)*[g(1,2) g(1,1)] ;
end

for i = 1 : rep2
    bn2(2,(i-1)*2+[1 2]) = (-1)^(i+ndy)*[g(2,2) g(2,1)] ;
end

for i = 1 : rep3
    bn3(2,(i-1)*2+[1 2]) = (-1)^(i)*[g(3,1) g(3,2)] ;
end



dl = h/ndy;
% Compute tension in bc 4
bc = zeros(1,2*ndy);
for i= 1:ndy
    pt1 = (ndx+1)*i;
    pt2 = (ndx+1)*i + ndx + 1;

    pt = ( cr(pt1,2)+cr(pt2,2) )/2;

    fx = -2*g(4,1)/h*pt + g(4,1);
    fy = 0.0;

    bc(1,2*(i-1)+1) = fx;
    bc(1,2*(i-1)+2) = fy;
end

for i = 1 : rep4
    bn4(2,(i-1)*2+[1 2]) = (-1)^(i+ndx)*[bc(1,2*(i-1)+1), bc(1,2*(i-1)+2)] ;
end


bn1(1,:) = bc1 ;
bn2(1,:) = bc2 ;
bn3(1,:) = bc3 ;
bn4(1,:) = bc4 ;

% Fix first element
bn3 = bn3(:,4:2:end);

%bn1 = [bc1; repmat([g(1,2),g(1,1)],1, rep1)] ;
%bn2 = [bc2; repmat([g(2,2),g(2,1)],1, rep2)] ;
%bn3 = [bc3; repmat(g(3,:),1, rep3)] ;
%bn4 = [bc4; repmat(g(4,:),1, rep4)] ;

return
