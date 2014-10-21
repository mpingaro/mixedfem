% by Marco Pingaro

function [stress,spost,rot] = solve(K,F,bc,ns,nd,nr)

% ASSIGN BOUNDARY CONDITION (sigma * n = 0)
for i = 1: size(bc,2)
    K(bc(1,i),:) = 0 ;
    K(bc(1,i),bc(1,i)) = 1 ;
end

% ASSIGN BOUNDARY CONDITION (sigma * n = t)
%
% Da fare!

%% SOLUTION OF LINEAR SYSTEM
soluz = K\F ;

nt = ns+nd+nr ;

stress= soluz(1:ns) ;
spost = soluz(ns+1:ns+nd) ;
rot = soluz(ns+nd+1:nt) ;

return
