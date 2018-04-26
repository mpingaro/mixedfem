% by Marco Pingaro

function [coordinates,element,mc] = beam(length,heigth,ndx,ndy)

%% Coorinates
px = length/ndx;
py = heigth/ndy;
npoint = (ndx+1)*(ndy+1);
xcor = [0:px:length];
ycor = [0:py:heigth];
for i=1:ndx+1
    y(i,:) = ycor;
end
for j=1:ndy+1
    x(:,j) = xcor; 
end
coordinates = zeros(npoint,1);
coordinates(:,1) = reshape(x,npoint,1);
coordinates(:,2) = reshape(y,npoint,1);

nelem = ndx*ndy;
element = zeros(nelem,4);
mc = zeros(nelem,13);
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
        mc((i-1)*ndx+j,[1 2 3 4 5 6 7 8 13]) = [i1 i2 i3 i4 i5 i6 i7 i8 i9];
    end
end
nedge = ndx*(ndy+1)+ndy*(ndx+1);
inz = 2*nedge;
fin = inz+4*nelem;
mc(:,9) = [inz+1:4:fin-3];
mc(:,10) = [inz+2:4:fin-2];
mc(:,11) = [inz+3:4:fin-1];
mc(:,12) = [inz+4:4:fin];

end
