function [am,aci,J,zm,zs] = bspline_compute(x,z,w,xm,dxm,lc,bctype);
% [am,aci,J,zm,zs] = bspline_compute(x,z,w,xm,dxm,lc,bctype);
% fit 1-D data to spline
%
% Input
%    x, the Nx1 data locations
%    z, the Nx1 observations
%    w, the Nx1 observation weights (e.g., rms(true)/(rms(error) + rms(true))
%    xm, the Mx1 EQUALLY SPACED grid nodes
%        xm(1) ought to be smaller than smallest x and xm(M) ought to be larger 
%        modify the xm or the input x to satisfy such that inconsistent use of bc does not flare up
%    dxm, the 1x1 grid spacing of xm (sure, could compute, just pass it in) used to scale
%    lc, the 1x1 spline curvature penalty weight
%        lc=4 wipes out wavelengths less than 2*dxm (Nyquist)
%    bctype, boundary condition type is either
%        2: second derivative vanishes at boundaries
%        1: first derivative vanishes at boundaries
%        0: value is zero at boundary
%        not specified: bc not enforced- good luck!
%
% Output
%    am, the Mx1 spline amplitudes
%    aci, the Mx1 error estimates
%    J, the 1x1 rms error of the spline surface
%    zm, the Mx1 spline values at the location xm
%    zs, the Nx1 spline values at the locations x
%      Note that spline values can be computed anywhere in the xm domain via
%       zs = bspline_curve(xi, xm, am, dxm);

% NOTE: if we are passing in gridded data to fix up, this code can be sped up by using the input grid indexes!!!!!

% input
N = length(x);
if(~exist('bctype'))
    bctype = -9;
end

% output grid
% d2 constraint performs well if boundary points entirely outside data
M = length(xm);

% normalize inputs
x = x/dxm;
xm = xm/dxm;
z = (w.*z);

% bc coeff
bc = sparse(M,1);
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
        bc(1:2) = [-4;1]; bc((M-1):M) = [-1;-4];
    otherwise
        % none: ;
end
% initial boundary (-1)
fb1 = bspline_basis(x-(xm(1)-1));
% end boundary (+1)
fbM = bspline_basis(x-(xm(M)+1));

% compute matrix coefficients
b = sparse(M,1); % data correlation
p = sparse(M,M); % model-model correlation at data
for m=1:M
    % forward sweep
    % set boundary function
    if(m<3)
        % initial boundary
        fb = fb1;
    elseif(m>=(M-3))
        % end boundary
        fb = fbM;
    else
        % interior
        fb = 0;
    end
    
    % supperpose boundary and central elements (or just central)
    if(m>1)
        % if we are past boundary, get the lagged terms that we have computed already
        f = g1;
    else
        % generate function
        f = bspline_basis(x-xm(m)) + bc(m)*fb; % values at data points    
    end
    b(m) = f'*(z); % spline-data covariance
    p(m,m) = f'*(w.*f); % spline-spline covariance, diagonal term
    
    % do first off diagonal terms
    if (m<M)
        mm = m+1;
        if(m>1 & m<(M-2))
            g1 = g2;
        else
            % include boundary influence here
            g1 = bspline_basis(x-xm(mm)) + bc(mm)*fb;
        end
        p(m,mm) = f'*(w.*g1);
        p(mm,m) = p(m,mm);
    end
    % do second off diagonal terms
    if (m<(M-1))
        mm = mm+1;
        if(m>1 & m<(M-2))
            g2 = g3;
        else
            g2 = bspline_basis(x-xm(mm)) +  bc(mm)*fb;
        end
        p(m,mm) = f'*(w.*g2);
        p(mm,m) = p(m,mm);
    end
    
   % do third off diagonal terms
    if (m<(M-2))
        mm = mm+1;
        g3 = bspline_basis(x-xm(mm)) +  bc(mm)*fb;
        p(m,mm) = f'*(w.*g3);
        p(mm,m) = p(m,mm);
    end
end
% q is the 2nd derivative covariance matrix
% it does not depend on the data
q = sparse(diag(ones(M,1))*6 + diag(ones(M-1,1),1)*(-27/8) + diag(ones(M-1,1),-1)*(-27/8) + diag(ones(M-3,1),3)*(3/8) + diag(ones(M-3,1),-3)*(3/8)); 

% implement the appropriate boundary conditions
q(1,1) = 3/4*bc(1)^2-9/4*bc(1)+3;
q(M,M) = q(1,1); % take advantage of symmetry

q(1,2) = -9/4+3/4*bc(1)*b(2)-9/8*bc(2);
q(2,1) = q(1,2);
q(M,M-1) = q(1,2);
q(M-1,M) = q(1,2);

q(1,3) = 3/8*bc(1);
q(3,1) = q(1,3);
q(M,M-2) = q(1,3);
q(M-2,M) = q(1,3);

q(2,2) = 3/4*bc(2)^2+21/4;
q(M-1,M-1) = q(2,2);

q(2,3) = -27/8+3/8*bc(2);
q(3,2) = q(2,3);
q(M-1,M-2) = q(2,3);
q(M-2,M-1) = q(2,3);

% compute the curvature penalty from lc
alpha = (lc/(2*pi))^4;

% this normalization allows q to kick in when sumw is small, else q not needed
if(length(w)==1)
    sumw = double(N*w + 1);
else
    sumw = double(sum(w)+1);
end

% add curvature terms
% this form enforces constant scaling of freq cutoff at about 0.25/lc
% try
    r = p/sumw + alpha*q/M;
% catch
%     lasterr
%     keyboard
% end
b = b/sumw;

% check matrix conditioning
% if no good, we can try for direct solution, but no error estimate
mrc = rcond(full(r));
if(mrc>100*eps)
%    fprintf('doing direct inversion...')
    r = inv(r);
    am = r*b;
else
%    fprintf('poor matrix condition, solving directly for am...')
    am = r\b;
    r = ones(size(r));
end
%dof = max([N-M, 1]);
%msz = (z'*z)/sumw;
%msm = (am')*(p)*(am)/M
%J = abs(msz-msm)
%aci = real(sqrt(diag(r)*J));

% here is code in regr_xzw
% msz = (z(id)'*z(id))/n;
% msr = msz - (b')*XX*(b);
% brmse = sqrt(diag(XX_inv)*msr/(n-m));

msz = double((z'*z)/sumw);
J = msz - (am')*(p/sumw + alpha*q/M)*am;
aci = sqrt(diag(r)*J/sumw);

% reconstruct
if(nargout>3)
    % evaluate at grid locations
    zm = bspline_curve(xm, xm, am, 1, bctype);
end
if(nargout>4)
% evaluate at data locations
zs = bspline_curve(x, xm, am, 1, bctype);
end
% done
