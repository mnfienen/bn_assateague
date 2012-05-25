function [x, y] = xyRotate(x,y,theta,xo,yo);

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
a = [cos(theta), sin(theta);-sin(theta), cos(theta)];

% translate to origin and rotate
x = [[x-xo],[y-yo]]*a;
y = x(:,2);
x = x(:,1);
