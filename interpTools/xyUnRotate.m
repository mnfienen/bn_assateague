function [x, y] = xyUnRotate(x,y,theta,xo,yo);

% rotate some data
%
% [xr, yr] = xyRotate(x,y,theta,xo,yo);
%
% Input
%   x,y are the un-rotated data 
%   theta is angle to rotate
%		(decimal degrees!)
%   xo,yo is the rotation origin
% 
% Output
%    xr, yr are rotated coords
%

% get rotation matrix
theta = pi*theta/180;
% unrotate
theta = - theta;

a = [cos(theta), sin(theta);-sin(theta), cos(theta)];



% translate to origin and rotate
x = [x,y]*a;
y = x(:,2)+yo;
x = x(:,1)+xo;
