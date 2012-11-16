function Zi = bspline_pertgrid(Zi,w, splinebctype, lc, dxm, dxi)
% [Zbsn] = bspline_pertgrid(Zi, w, splinebctype, lc, dxm, dxi)
% 
% return a cubic-b-splined version of input Zi using weights, w
% assume that we are splining perturbations, that are driven to zero at boundaries
% the method is sensitive to the weights.  If weights are zero, solution heads for bc 
%
% Input
%   Zi, the data (assume Nx1 or NxM regular gridded data)
%   w, weights for the data
%   e.g. [wt, var_z] = consistentWeight(z, s.^2, wtol); 
%      OR w = (1-NEi)+delta_e;
%         where delta_e is some small number like 10-9, which indicates how strongly to drive the solution to zero where no data exist
%   splinebcytype
%        2: second derivative vanishes at boundaries
%        1: first derivative vanishes at boundaries
%        0: value is zero at boundary
%       10: force to value and derivative to zero
%        not specified: bc not enforced- good luck!
%   lc, the smoothing constraint value (usually 4, to eliminate aliasing)
%   dxm, the coarsening of the grid (e.g., 2 means compute with dx that is twice input dx)
%   dxi, fining (e.g., 0.1 means return spline on grid that is ten times finer than input dx)
%
% Output
%  Zbsn, the bsplined grid
[Ny,Nx,Nd] = size(Zi);
% allow 1-d input
if(Ny>1 & Nx==1)
    Zi = Zi';
    w = w';
    [Ny,Nx,Nd] = size(Zi);
end
id = find(isnan(w));
w(id) = 0;

% set penalty function for curvature
if(~exist('lc')==1)
    lc=4;
end

% specify bc condition
fix_ends = 0;
if(~exist('splinebctype')==1)
    splinebctype = zeros(1,Nd); % value zero boundary
elseif(Nx>1 & Ny>1 & length(splinebctype)==1)
    % make bc type variable with each dimension
    splinebctype(2) = splinebctype(1);
    fix_ends = [0,0];
elseif(length(splinebctype)==2)
    fix_ends = [0,0];
end

% and check to see if we are pinning bc
for i=1:length(fix_ends)
    if(splinebctype(i)==10)
        splinebctype(i) = 1; % compute as if using derivative constraint
        fix_ends(i) = 1;
    end
end

% specify coarseness of compute grid
if(~exist('dxm')==1)
    dxm = 2; % make it coarser than input
end

% input grid
x = 1:Nx;

% outputgrid
if(~exist('dxi')==1)
    dxi = 1;
end
xi = 1:dxi(1):Nx;
Nxi=length(xi);

% put xm on boundary
nxm = fix([x(1,end)-x(1,1)]/dxm(1));
dxm(1) = [x(1,end)-x(1,1)]/nxm;
xm = x(1,1)+dxm(1)*[0:nxm]';

% can proceed in 1-d, or...
% check for 2-d
if(Ny>1 & Nx>1)
    y = [1:Ny]';
    % repeat for all dims
    %if(length(dxm>1))
    if(length(dxm)>1)
        dym = dxm(2);
        dxm = dxm(1);
    else
        dym = dxm;
    end
    nym = fix([y(end,1)-y(1,1)]/dym);
    if(nym<3)
        nym = 3; % 29 dec 2010: implement this to support 2-d anaylyiss for all cases
    end
    dym = [y(end,1)-y(1,1)]/nym;
    ym = y(1,1)+dym*[0:nym]';
    if(length(dxi)>1)
        dyi = dxi(2);
        dxi = dxi(1);
    else 
        dyi = dxi;
    end
    yi = [1:dyi:Ny]';
    Nyi = length(yi);
    
    Am = zeros(length(ym),length(xm));
    Cm = zeros(length(ym),length(xm));
    % spline alongshore first
    ztmp0 = zeros(size(Zi(:,1)));
    for i=1:Nx
        fprintf('.')
        id = find(~isnan(Zi(:,i)) & w(:,i)>eps);
        ztmp = ztmp0;
        ztmp(id) = Zi(id,i);
        try
        [am,aci,J] = bspline_compute(y,ztmp,w(:,i),ym,dym,lc, splinebctype(2));
        catch
            keyboard
        end
        Am(:,i) = am;
        aci(aci==J)=nan;
        Cm(:,i) = aci/max(aci);
    end
    if(fix_ends(2))
        Am(1,:)=-Am(2,:)/2;
        Am(end,:)=-Am(end-1,:)/2;
        Cm(1,:)=0;
        Cm(end,:)=0;    
    end
    Zi = Am; % spline the Ams
    % update weights
    minw = min(w(:));
    w  = minw + 1-Cm./(1+Cm);
    w(isnan(w))=(1/length(w));
    yprime = y; % save this
    y = ym; % now compute on this grid
    Ny = length(ym);
end
fprintf('\n')

% now run the spline along the x direction
Zprime = zeros(Ny,Nxi);
ztmp0 = zeros(size(Zi(1,:)));
for i=1:Ny
    % run spline
    fprintf('o')
    id = find(~isnan(Zi(i,:)) & w(i,:)>eps);
    ztmp = ztmp0;
    ztmp(id) = Zi(i,id);
    [am] = bspline_compute(x',ztmp',w(i,:)',xm,dxm,lc, splinebctype(1));
    % check if we are also forcing zero at boundary
    if(fix_ends(1))
        am(1)=-am(2)/2;
        am(end)=-am(end-1)/2;
    end
    % get result on fine grid
    zm = bspline_curve(xi', xm, am, dxm, splinebctype(1))';
    Zprime(i,:) = zm;
end
fprintf('\n')
Zi = Zprime;

% now convert back, Zi are really Am^y
if(Ny>1 & Nx>1)
    y = yprime;
    Ny = length(y);
    Ziprime = Zi;
    % check if we are also forcing zero at boundary
    if(fix_ends)
        Ziprime(1,:)=-Ziprime(2,:)/2;
        Ziprime(end,:)=-Ziprime(end-1,:)/2;
    end
    Zi = zeros(Nyi,Nxi);
    for i=1:Nxi
        fprintf('+')
        zm = bspline_curve(yi, ym, Ziprime(:,i), dym, splinebctype(2));
        Zi(:,i) = zm;
    end
end
fprintf('\n')
