% by Marco Pingaro

function plotsol(coordinates,element,stress,spost,rot,def,ndx,ndy)

% PLOT SOLUTION DISPLACEMENT   
x = reshape(coordinates(:,1),ndx+1,ndy+1) ;
y = reshape(coordinates(:,2),ndx+1,ndy+1) ;
def_x = reshape(def(:,1),ndx+1,ndy+1) ;
def_y = reshape(def(:,2),ndx+1,ndy+1) ;

spost_x = spost(1:6:6*ndx*ndy-5) ;
spost_y = spost(4:6:6*ndx*ndy) ;
rot = reshape(rot,ndx+1,ndy+1);

% Undeformed mesh
figure, surf(x,y,zeros(ndx+1,ndy+1))
axis equal
view(0,90)

% Deformed mesh
figure, surf(def_x,def_y,zeros(ndx+1,ndy+1))
axis equal
view(0,90)

% Displacemnt X
figure, surf(x,y,def_x)
axis equal
view(60,80)

% Displacemnt Y
figure, surf(x,y,def_y)
axis equal
view(60,80)

% Displacement X
figure
for j = 1:ndx*ndy
    cx = coordinates(element(j,1:4),1) ;
    cy = coordinates(element(j,1:4),2) ;
    sp_x = spost_x(j)*ones(4,4) ;
    surf(cx,cy,sp_x)
    hold on
end
axis equal
view(60,80)

% Displacement Y
figure
for j = 1:ndx*ndy
    cx = coordinates(element(j,1:4),1) ;
    cy = coordinates(element(j,1:4),2) ;
    sp_y = spost_y(j)*ones(4,4) ;
    surf(cx,cy,sp_y)
    hold on
end
axis equal
view(60,80)

% Rotation
figure, surf(x,y,rot)
axis equal
view(60,80)

% Stress 
% da implementare

return
