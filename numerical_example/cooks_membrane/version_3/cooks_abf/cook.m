% by Marco Pingaro

function [coordinates,element,mc] = cook(nodes,ndx,ndy,dl1,dl2)

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

nelem = ndx*ndy;
element = zeros(nelem,4);
mc = zeros(nelem,15);
for i =1:ndy
    for j = 1:ndx
        %% Element
        element(ndx*(i-1)+j,1) = (ndx+1)*(i-1)+j;
        element(ndx*(i-1)+j,2) = element(ndx*(i-1)+j,1)+1;
        element(ndx*(i-1)+j,3) = (ndx+1)*i+j+1;
        element(ndx*(i-1)+j,4) = element(ndx*(i-1)+j,3)-1;
        %% Corrispondence Matrix
        i1 = (4*ndx+2)*(i-1)+2*(j-1)+1;
        i2 = i1+1;
        i7 = (4*ndx+2)*(i-1)+2*(j-1)+1+2*ndx;
        i8 = i7+1;
        i3 = i8+1;
        i4 = i3+1;
        i5 = (4*ndx+2)*(i-1)+4*ndx+2*(j-1)+3;
        i6 = i5+1;
        i9 = (-1)^(i+j);
        mc((i-1)*ndx+j,[1 2 3 4 5 6 7 8 15]) = [i1 i2 i3 i4 i5 i6 i7 i8 i9];
    end
end
nedge = ndx*(ndy+1)+ndy*(ndx+1);
inz = 2*nedge;
fin = inz+6*nelem;
mc(:,9) = [inz+1:6:fin-5];
mc(:,10) = [inz+2:6:fin-4];
mc(:,11) = [inz+3:6:fin-3];
mc(:,12) = [inz+4:6:fin-2];
mc(:,13) = [inz+5:6:fin-1];
mc(:,14) = [inz+6:6:fin];

return
