%%
nodes   = [0, 0; 48, 44; 48, 60; 0, 44] ;
dl1     = nodes(3,2)-nodes(2,2) ;
dl2     = nodes(4,2) ;

ndx = 64;
ndy = 64;

[coordinates,element] = GenMeshCook(nodes,ndx,ndy,dl1,dl2);
plotsol(coordinates,ndx,ndy);


function [coordinates,element] = GenMeshCook(nodes,ndx,ndy,dl1,dl2)

%% Coorinates
% theta
theta = zeros(ndy+1,1);
for i = 0:length(theta)-1
    theta(i+1,1) = atand((nodes(2,2)+i*dl1/ndy-i*dl2/ndy)/nodes(2,1));
end

%% Coordinates Matrix
% This matrix contains at rows the nodal points of physical element and at 
% columns the relative coordinates
coordinates = zeros((ndx+1)*(ndy+1),2);
for i = 1:ndy+1
    for j = 1:ndx+1
        coordinates(j+(i-1)*(ndx+1),1) = (j-1)*nodes(2,1)/ndx;
        coordinates(j+(i-1)*(ndx+1),2) = coordinates(j+(i-1)*(ndx+1),1)*tand(theta(i,1))...
        +(i-1)*nodes(4,2)/ndy;
    end
end

var = zeros((ndx+1)*(ndy+1),2);
for i=1:ndy-1
    for j=1:ndx-1
        ind =(ndx+1)*(i-1)+ndx+2+j;
        var(ind,1) = (-1)^(i*j)*rand(1,1)*coordinates(ind,1)/coordinates(ind+1,1)/10;
        var(ind,2) = (-1)^(i*j)*rand(1,1)*4*coordinates(ind,2)/coordinates(ind+1,2)/20;
    end
end
coordinates = coordinates+var;


nelem = ndx*ndy;
element = zeros(nelem,4);
for i =1:ndy
    for j = 1:ndx
        %% Element
        element(ndx*(i-1)+j,1) = (ndx+1)*(i-1)+j;
        element(ndx*(i-1)+j,2) = element(ndx*(i-1)+j,1)+1;
        element(ndx*(i-1)+j,3) = (ndx+1)*i+j+1;
        element(ndx*(i-1)+j,4) = element(ndx*(i-1)+j,3)-1;
    end
end

end


function plotsol(coordinates,ndx,ndy)

% PLOT SOLUTION DISPLACEMENT   
x = reshape(coordinates(:,1),ndx+1,ndy+1) ;
y = reshape(coordinates(:,2),ndx+1,ndy+1) ;

% Undeformed mesh
figure, surf(x,y,zeros(ndx+1,ndy+1))
axis equal
view(0,90)

end