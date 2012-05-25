function [z, e] = bspline_curve(y, ym, am, dx, bctype, aci);
% [zi, ei] = bspline_curve(xi, xm, am, dxm, bctype aci);
% generate the spline surface from a series of evaluated basis functions
%
% Input
%   xi, the Nx1 locations of interst
%   xm, the Mx1 EQUALLY SPACED locations of the basis function
%   am, the corresponding basis function amplitudes
%   dxm, the grid spacing of xm
%     bctype, boundary condition type is either
%         2: second derivative vanishes at boundaries
%         1: first derivative vanishes at boundaries
%         0: value is zero at boundary
%         not specified: bc not enforced- good luck!%
% Output
%    zi, the Nx1 spline values at each xi
%
% Notes
% to enforce both zero at boundary and zero slope, choose bctype=1 and set am(1)=-0.5*am(2); 

M = length(ym);
if(exist('bctype'))
    % bc coeff
    bc = zeros(M,1);
    switch bctype
        % bc effective if consistent with the data
        case 2
            % d2=>0: 
            bc(1:2) = [2;-1]; bc((M-1):M) = [-1;2];
        case 1
            % d1=>0: 
            bc(1:2) = [0;1]; bc((M-1):M) = [1;0];
        case 0
            % d0=>0: 
            bc(1:2) = [-4;-1]; bc((M-1):M) = [-1;-4];
        otherwise
            % none: ;
    end
end

% normalize by length scale, dx
y=y/dx;
ym = ym/dx;
z = zeros(size(y));
e = zeros(size(y));
% compute basis function value at each output location and mult by amplitude
for i=1:length(y)
    % use only those functions that support location y(i)
    % saves computing alot of zeros
    id = find(abs(y(i)-ym)<=2);
    if(length(id)>0)
        f = bspline_basis(y(i)-ym(id));
        z(i) = am(id)'*f;
        if(exist('aci'))
            e(i) = ((aci(id)').^2)*(f.^2);
        end
        if(any(id<3) & exist('bctype'))
            % include bc(1)
            abc = am(1)*bc(1) + am(2)*bc(2);
            f1 = bspline_basis(y(i)-ym(1)+1);
            z(i) = z(i) + abc*f1; 
            if(exist('aci'))
                e(i) = e(i) + ((aci(1)*bc(1)).^2 + (aci(2)*bc(2)).^2)*(f1^2);
            end
        end
        if(any(id>(M-2)) & exist('bctype'))
            % include bc(M)
            abc = am(M-1)*bc(M-1) + am(M)*bc(M);
            f2 = bspline_basis(y(i)-ym(M)-1);
            z(i) = z(i) + abc*f2; 
            if(exist('aci'))
                e(i) = e(i) + ((aci(M-1)*bc(M-1)).^2 + (aci(M)*bc(M)).^2)*(f2^2);
            end
        end
    else
        y(i) = 0;
    end
end