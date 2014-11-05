% by Marco Pingaro

function [spost,pres] = solve(K,F,bc,nd,np)

% ASSIGN BOUNDARY CONDITION (u = 0)
for i = 1: size(bc,2)
    K(bc(1,i),:) = 0 ;
    K(bc(1,i),bc(1,i)) = 1 ;
    F(bc(1,i),1) = 0 ;
end

% ASSIGN BOUNDARY CONDITION (p = t)
%
% Da fare!

%% SOLUTION OF LINEAR SYSTEM
soluz = K\F ;

nt = nd+np ;

spost = soluz(1:nd) ;
pres  = soluz(nd+1:nt) ;

return
