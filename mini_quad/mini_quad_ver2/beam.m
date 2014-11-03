% by Marco Pingaro

function [coordinates,element,bc1,bc2,bc3,bc4] = beam(length,heigth,ndx,ndy)

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
for i =1:ndy
    for j = 1:ndx
        %% Element
        element(ndx*(i-1)+j,1) = (ndx+1)*(i-1)+j;
        element(ndx*(i-1)+j,2) = element(ndx*(i-1)+j,1)+1;
        element(ndx*(i-1)+j,3) = (ndx+1)*i+j+1;
        element(ndx*(i-1)+j,4) = element(ndx*(i-1)+j,3)-1;
    end
end

bl1 = 2.*[1:ndx+1] ;
bl2 = 2.*[ndy*(ndx+1)+1:npoint] ;
bl3 = 2.*[1:ndx+1:(ndx+1)*ndy+1] ;
bl4 = 2.*[ndx+1:ndx+1:npoint] ;

bc1 = [bl1-1, bl1];
bc2 = [bl2-1, bl2];
bc3 = [bl3-1, bl3];
bc4 = [bl4-1, bl4];

return
