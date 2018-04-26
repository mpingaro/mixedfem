% Written by Marco Pingaro and Paolo Venini
%

function [W,P] = GaussQuad2D(degx,degy)

    % Call simple 1D gauss quadrature in the two directions
    [Px, Wx] = GaussQuad(degx);
    [Py, Wy] = GaussQuad(degy);

    % Number of points per directions
    nx = size(Px,1);
    ny = size(Py,1);

    %%
    P = zeros(nx*ny, 2);
    W = zeros(nx*ny, 1);

    k=1;
    for i=1:nx
        for j=1:ny
            W(k,1) = Wx(i)*Wy(j);
            P(k,[1, 2]) = [Px(i), Py(j)];
            k = k + 1;
        end
    end

return
