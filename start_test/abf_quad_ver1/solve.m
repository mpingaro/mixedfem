% by Marco Pingaro

function [stress,spost,rot] = solve(K,F,bn,g,ndx,ndy,ns,nd,nr)

% ASSIGN BOUNDARY CONDITION (sigma * n = t)
for iside = bn
    cd = neumann(ndx,ndy,g(iside,:),iside);
    for ndof = 1: size(cd,2)
        K(cd(1,ndof),:) = 0 ;
        K(cd(1,ndof),cd(1,ndof)) = 1 ;
        F(cd(1,ndof),1) = cd(2,ndof) ;
    end
end

%% SOLUTION OF LINEAR SYSTEM
soluz = K\F ;

nt = ns+nd+nr ;

stress= soluz(1:ns) ;
spost = soluz(ns+1:ns+nd) ;
rot = soluz(ns+nd+1:nt) ;

return
