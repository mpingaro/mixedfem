% by Marco Pingaro

function [stress,spost,rot] = solve(K,F,bn,ns,nd,nr)

% ASSIGN BOUNDARY CONDITION (sigma * n = t)
% for i = 1:size(bn,2)
%     K(bn(1,i),:) = 0 ;
%     K(bn(1,i),bn(1,i)) = 1 ;
%     F(bn(1,i),1) = bn(2,i) ;
% end

%% SOLUTION OF LINEAR SYSTEM
soluz = K\F ;

nt = ns+nd+nr ;

stress= soluz(1:ns) ;
spost = soluz(ns+1:ns+nd) ;
rot = soluz(ns+nd+1:nt) ;

return
