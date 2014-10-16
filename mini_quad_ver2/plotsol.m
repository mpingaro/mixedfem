% by Marco Pingaro
% PLOT SOLUTION 

function plotsol(coordinates,spost,pres,ndx,ndy,nnod)

x = reshape(coordinates(:,1),ndx+1,ndy+1) ;
y = reshape(coordinates(:,2),ndx+1,ndy+1) ;
spost_x = spost(1:2:2*nnod-1) ;
spost_y = spost(2:2:2*nnod) ; 
spost_x = reshape(spost_x,ndx+1,ndy+1) ;
spost_y = reshape(spost_y,ndx+1,ndy+1) ;
def_x = x+spost_x ;
def_y = y+spost_y ;
pres = reshape(pres,ndx+1,ndy+1) ;

% Undeformed mesh
figure, surf(x,y,zeros(ndx+1,ndy+1))
axis equal
view(0,90)

% Deformed mesh
figure, surf(def_x,def_y,zeros(ndx+1,ndy+1))
axis equal
view(0,90)

% X directions
figure, surf(x,y,spost_x)
hold on
axis equal
view(60,80)

% Y directions
figure, surf(x,y,spost_y)
hold on
axis equal
view(60,80)

% Pression
figure, surf(x,y,pres)
hold on
axis equal
view(60,80)

end
