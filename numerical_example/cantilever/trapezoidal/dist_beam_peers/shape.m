%% --------------------------------------------------------------------- %%
% ----------------------------------------------------------------------- % 
% --- SHAPE FUNCTION RT_0 QUAD ( by Marco Pingaro & Giacomo Maurelli      %
%                                                                         %
%                                                                         %
%                                                                         %
% Stress       : RT_0 = [a + b*x ; c + d*y]                 Continuous    %
% Stress       : Stabilization [x^2 -1 ; 0] , [0 ; y^2 - 1] Discontinuous %
% Displacement : P_0 = [a ; b]                              Discontinuous %
% Rotation     : P_1 bi-linear = a + b*x + c*y + d*xy       Continuous    %
%                                                                         %
% Element Domain (-1, 1) x (-1, 1)                                        %
%                                                                         %
% ----------------------------------------------------------------------- %
%% SHAPE

%% -- Stress
%sigma(1) = [0; -0.5+0.5*y];
%sigma(2) = [0.5+0.5*x; 0];
%sigma(3) = [0; 0.5+0.5*y];
%sigma(4) = [-0.5+0.5*x; 0]; 
%sigma(5) = [x*x - 1; 0];
%sigma(6) = [0; y*y - 1];
%% --- Divergence
%div(1) = 0.5;
%div(2) = 0.5;
%div(3) = 0.5;
%div(4) = 0.5;
%div(5) = 0;
%div(6) = 0;
%% -- Displacement
%disp(1) = [1; 0];
%disp(2) = [0; 1];
%% -- Rotation
%rot(1) =  0.25*(1-x)*(1-y);
%rot(2) =  0.25*(1+x)*(1-y);
%rot(3) =  0.25*(1+x)*(1+y);
%rot(4) =  0.25*(1-x)*(1+y);
