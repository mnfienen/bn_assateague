function [am,aci,J,b,r] = bspline_compute(x,y,z,w,xm,ym,lc,bctype);
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
x = x(:);
y = y(:);
z = z(:);
w = w(:);

% output grid
% d2 constraint performs well if boundary points entirely outside data
Mx = length(xm);
My = length(ym);

% normalize inputs
dxm = diff(xm(1:2));
dxm = [dxm(1);diff(xm(:))];
if(My>1)
    dym = diff(ym(1:2));
    dym = [dym(1);diff(ym(:))];
else
    dym=1;
end
x = x./interp1(xm,dxm,x);
if(My>1)
    y = y./interp1(ym,dym,y);
else
    y = y/dym;
end
xm = xm./dxm;
ym = ym./dym;
z = (w.*z);

% pre-compute the curvature terms
if(lc~=0)
    fprintf('precomputing curvature correlations...\n')
    dxtmp = 0.1;
    xtmp = -4:dxtmp:4;
    mtmp = length(xtmp);
    [xtmp,ytmp] = meshgrid(xtmp,xtmp);
    f = reshape(bspline_basis(xtmp(:)).*bspline_basis(ytmp(:)),size(xtmp));
    [gx,gy] = gradient(f); % don't include the stepsize here, so it is assumed to be one, so divide by stepsize
    [gxx,gxy] = gradient(gx); % divide by stepsize again for 2nd derivatative, so missing dx^2 factor
    [gyx,gyy] = gradient(gy);
    % convolution is integral, so mult by step size (dx*dy), this cancels the missing stepsizes!
    Q= conv2(gxx,gxx,'same') + conv2(gyy,gyy,'same') + conv2(gxy,gxy,'same') + conv2(gyx,gyx,'same');
    % get elements at spline nodes
    qp = Q(1:(1/dxtmp):end, 1:(1/dxtmp):end);
    clear Q gx gy gxx gyy gxy gyx junk xtmp ytmp mtmp dxtmp f 
else
    qp = zeros(9,9);
end

% set bc constraints
bcx = sparse(Mx,1);
bcy = sparse(My,1);
% for now, dont mix and match
switch bctype
    % bc effective if consistent with the data
    case 2
        % d2=>0: 
        bcx([1;2;Mx-1;Mx]) = [2;-1;-1;2]; 
        bcy([1;2;My-1;My]) = [2;-1;-1;2]; 
    case 1
        % d1=>0: 
        bcx([1;2;Mx-1;Mx]) = [0;1;1;0]; 
        bcy([1;2;My-1;My]) = [0;1;1;0]; 
    case 0
        % d0=>0: 
        bcx([1;2;Mx-1;Mx]) = [-4;1;-1;-4];
        bcy([1;2;My-1;My]) = [-4;1;-1;-4];
    otherwise
        % none: ;
end

% compute matrix coefficients
fprintf('computing correlations...\n')
M = Mx*My;
b = sparse(M,1); % data correlation at xm,ym
p = sparse(M,M); % model-model correlation at xm,ym
q = p; % model-model 2nd derivative correlation at xm,ym
hw = waitbar(0,'bspline2d');
tstart = cputime;
for mx=1:Mx
    waitbar((mx-1)*My/M, hw, sprintf('bspline2d: %.3f seconds', cputime-tstart)); 
    % speed things up by selecting relevant data now (data can touch the boundary spline)
    xid = find(abs(x-xm(mx))<3);
    if(length(w)==1)
        % weights are same everywhere, so pass this index
        wxid = 1;
    else
        % weights given at each data point
        wxid = xid;
    end
    if(mx<=2)
        % set it now
        fbx = sparse(bspline_basis(x(xid)-(xm(1)-1))); % don't forget this has a y-dim!
    elseif(mx>=(Mx-1))
        % reset it
        fbx = sparse(bspline_basis(x(xid)-(xm(Mx)+1)));
    else
        fbx = 0;
    end
    fx = sparse(bspline_basis(x(xid)-xm(mx))) + bcx(mx)*fbx;
    for my=1:My
        % here is index into Mx*My-by-Mx*My matrix
        id = (mx-1)*My+my;
        % supperpose boundary and central elements (or just central)
        if(my<=2)
            % set it
            fby = sparse(bspline_basis(y(xid)-(ym(1)-1)));
        elseif(my>=(My-1))
            % reset it
            fby = sparse(bspline_basis(y(xid)-(ym(My)+1)));
        else
            fby = 0;
        end
        fy = sparse(bspline_basis(y(xid)-ym(my))) + bcy(my)*fby;
        f = fx.*fy;
        b(id) = f'*(z(xid)); % spline-data covariance
        
        % get covariances
        for ixprime = 0:3
            if(ixprime==0)
                % fx unchanged
                fxprime = fx;
            elseif((mx+ixprime)<=Mx) % ixprime is 1-sided
                %get the correct boundary function at ALL lags
                fxprime = sparse(bspline_basis(x(xid)-xm(mx)-ixprime)) + bcx(mx+ixprime)*fbx;
            else
                fxprime = 0;
            end
            for iyprime = -3:3
                % lagged index
                idprime = (mx+ixprime-1)*My+(my+iyprime);
                % now, check here that we are on domain, or do nothing if done already cause I can't figure out indices
                if((mx+ixprime)>=1 & (mx+ixprime)<=Mx & (my+iyprime)>=1 & (my+iyprime)<=My & p(id, idprime)==0) 
                    if(iyprime==0)
                        % fy unchanged
                        fyprime = fy;
                    elseif((my+iyprime)>=1 & (my+iyprime)<=My) % iyprime is 2-sided
                        fyprime = sparse(bspline_basis(y(xid)-ym(my)-iyprime)) + bcy(my+iyprime)*fby;
                    else
                        fyprime = 0;
                    end
                    % spline-spline covariance
                    fprime = fxprime.*fyprime;
                    pprime = f'*(w(wxid).*fprime);
                    p(id, idprime) = pprime;
                    p(idprime, id) = pprime;
                    % curvature-curvature covariance
                    qprime = qp(5+ixprime,5+iyprime);
                    q(id, idprime) = qprime; 
                    q(idprime, id) = qprime; 
                end
            end
        end
    end
end
tend = cputime;
fprintf('elapsed time is %.3f seconds\n', tend-tstart)
delete(hw)

% compute the curvature penalty from lc
alpha = (lc/(2*pi))^4;

% this normalization allows q to kick in when sumw is small, else q not needed
if(length(w)==1)
    sumw = N*w + 1;
else
    sumw = sum(w)+1;
end
b = b/sumw;

% add curvature terms
% this form enforces constant scaling of freq cutoff at about 0.25/lc
r = p/sumw + alpha*q/(Mx*My);

% check matrix conditioning
try 
    MYMXMAX = 1000;
    if(My*Mx>MYMXMAX)
        % we can trigger iterative methods this way
        error(sprintf('array exceeds max dimension %d ...', MYMXMAX))
    end
    fprintf('trying direct inversion\n')
    warning('');
    am = inv(r)*(b);
    aci = inv(r)*ones(My*Mx,1);
    if(length(findstr('singular', lastwarn))>0)
        % poor conditioning, try this way 
        fprintf('slash inversion\n')
        am = r\(b);
        aci = r\ones(My*Mx,1);
    end
catch
    % fails only if out of memory
    disp(lasterr)
    fprintf('trying iterative inversion\n')
    tol = [];
    maxit = 500;
    am = cgs(r,(b),tol,maxit); %works, fast too
    aci = cgs(r,ones(Mx*My,1),tol, maxit);
end

% scale error
msz = (z'*z)/sumw;
J = msz - (am')*(p/sumw + alpha*q/(Mx*My))*am;
aci = sqrt(aci*J/sumw);

% fixup output
am = reshape(am, My,Mx);
aci = real(reshape(aci, My,Mx));
