% by Marco Pingaro

function [gauss_p,gauss_w,npg] = quadrature_16()

%% Gauss Quadrature (16 points)
%
qp = [-0.8611363115940525, ...
    -0.3399810435848563, 0.8611363115940525, 0.3399810435848563];

wp = [0.3478548451374544, ... 
    0.6521451548625460, 0.3478548451374544, 0.6521451548625460];

npg = 16; 
gauss_p = zeros(2,npg);
gauss_w = zeros(1,npg);
for i = 1:4
    for j = 1:4
        gauss_p(1, 4*(i-1)+j) = qp(i);
        gauss_p(2, 4*(i-1)+j) = qp(j);
        gauss_w(1, 4*(i-1)+j) = wp(i)*wp(j);
    end
end

return
