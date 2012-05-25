function [Zbsn] = bspline_grid(x,y,Zi,w)
% [Zbsn] = bspline_grid(x,y,Zi,w)
% return a cubic-b-splined version of input Zi using weights, w
%   x,y, NxM grid of coordinate values
%   Zi, the data
%   w, weights for the data
%   e.g. [wt, var_z] = consistentWeight(z, s.^2, wtol); OR w = (1-NEi)+eps;
% corrected some errors 

[Ny,Nx,Nd] = size(x);
id = find(isnan(w));
w(id) = 0;

% init
Zbsn = repmat(nan, Ny,Nx);
zmin = min(-Zi(:));
zmax = max(-Zi(:));


% set penalty function for curvature
lc=4;

% specify bc condition
splinebctype = 1; % first derivative zero at boundary

% output grid
% d2 constraint performs well if boundary points entirely outside data
dx = diff(x(1,1:2));
dxm = 2*dx; % make it coarser than input
% get points
xm = [x(1,1):dxm:x(1,end)]';
xm = [xm(1)-dxm; xm; xm(end)+dxm];

% can proceed in 1-d, or...

% check for 2-d
if(Ny>1 & Nx>1)
    % repeat for all dims
    dy = diff(y(1:2,1));
    dym = 2*dy;
    ym = [y(1,1):dym:y(end,1)]';
    ym = [ym(1)-dym; ym; ym(end)+dym];
    
    % spline each cross-section
    % spline alongshore first
    for i=1:length(x(1,:))
        fprintf('spline %d of %d\n', i, Nx)
        % run spline
        id = find(~isnan(Zi(:,i))  & y(:,i)>y(1,i) & y(:,i)<y(end,i));
        ztmp = Zi(id,i);
        slope = (ztmp(end) - ztmp(1))/(y(id(end),i) - y(id(1),i));
        ztrend = slope* y(:,i) - slope*y(id(1),i) + Zi(id(1),i); 
        ztmp = ztmp-ztrend(id);
        [am,aci,J] = bspline_compute(y(id,i),ztmp,w(id,i),ym,dym,lc, splinebctype);
        zm = bspline_curve(y(:,1), ym, am, dym);
        Zbsn(:,i) = zm + ztrend;
    end
    % put output into input
    Zi = Zbsn;
    
end

% now run the spline along the x direction
Zbsn = repmat(nan, Ny,Nx);
for i=1:length(y(:,1))
    % run spline
    fprintf('spline %d of %d\n', i, Ny)
    id = find(~isnan(Zi(i,:)) & x(i,:)>x(i,1) & x(i,:)<x(i,end));
    ztmp = Zi(i,id);
    slope = (ztmp(end) - ztmp(1))/(x(i,id(end)) - x(i,id(1)));
    ztrend = slope* x(i,:) - slope*x(i,id(1)) + ztmp(1); 
    ztmp = ztmp-ztrend(id);
    [am,aci,J] = bspline_compute(x(i,id)',ztmp',w(i,id)',xm,dxm,lc, splinebctype);
    zm = bspline_curve(x(1,:)', xm, am, dxm)';
    Zbsn(i,:) = zm + ztrend;
end

