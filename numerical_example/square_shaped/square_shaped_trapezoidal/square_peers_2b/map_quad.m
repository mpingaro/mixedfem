%% Map square element

% [b1, c1;b2, c2]*[xi, eta] = [d1 - a1*xi*eta, d2 - a2*xi*eta]
% d1 = 4*x - (x1+x2+x3+x4)
% d2 = 4*y - (y1+y2+y3+y4)
%
% coef = [a1, a2; b1, b2; c1, c2]
% Example
%
% Physic points
% pt = [0,0;
%     5,0;
%     5,5;
%     0,5];
% Point in reference domain
% ev_pt = [1,-1];

function [x,y] = map_quad(pt,ev_pt)

% Counterclockwise numeration of node
mat = [1, -1, 1, -1;
    -1, 1, 1, -1;
    -1, -1, 1, 1];

coef = mat*pt;
sumx = sum(pt(:,1));
sumy = sum(pt(:,2));

x = ( coef(2,1)*ev_pt(1)+coef(3,1)*ev_pt(2) + coef(1,1)*ev_pt(1)*ev_pt(2) + sumx )/4;
y = ( coef(2,2)*ev_pt(1)+coef(3,2)*ev_pt(2) + coef(1,2)*ev_pt(1)*ev_pt(2) + sumy )/4;

end