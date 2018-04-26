% Written by Marco Pingaro and Paolo Venini
%
clear; close all; clc;

syms x y l m;
% Displacement 
u(1) = cos(pi.*x)*sin(2*pi.*y);
u(2) = sin(pi.*x)*cos(pi.*y);

% Strain tensor
epsi(1,1) = diff(u(1), x);
epsi(2,2) = diff(u(2), y);
epsi(1,2) = 1/2.*( diff( u(1), y ) + diff( u(2), x ) );
epsi(2,1) = epsi(1,2);

% Divergence of displacement
divu = diff( u(1), x ) + diff( u(2), y ); 

% Stress tensor
sig = 2.*m.*epsi + l.*divu.*blkdiag(1,1);

% Compute the divergence
div(1) = diff( sig(1,1), x ) + diff( sig(1,2), y );
div(2) = diff( sig(2,1), x ) + diff( sig(2,2), y );

disp( div(1) );
disp( div(2) );
