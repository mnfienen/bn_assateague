function [y] = bspline_basis(y);
% [y] = bspline_basis(x);
% b-spline basis function evaluated at x
%
% Input
% x, the Nx1 independent variable locations
%     x is normalized by dx (grid node spacing) and centered on xm (node of interest)
%
% Output
%    y, the basis function for x-xm=0

% function is symmetric
ya = abs(y);
y=zeros(size(y));

% and defined piecewise
id12 = find(1<=ya & ya<2);
id1 = find(ya<1);
y(id12) = 0.25*((2-ya(id12)).^3);
y(id1) = 0.25*((2-ya(id1)).^3)-((1-ya(id1)).^3);
