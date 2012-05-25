function [z, e] = bspline_curve2d(x,y, xm,ym, am, bctype, aci);
% [zi, ei] = bspline_curve2d(x,y, xm,ym, am, bctype, aci);
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

My = length(ym);
Mx = length(xm);

% set bc constraints
bcx = sparse(Mx,1);
bcy = sparse(My,1);
% for now, dont mix and match
if(exist('bctype'))
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
end

% normalize by length scale, dx
dxm = diff(xm(1:2));
dxm = [dxm(1);diff(xm(:))];
if(My>1)
    dym = diff(ym(1:2));
    dym = [dym(1);diff(ym(:))];
else
    dym=1;
end
%x = x./interp1(xm,dxm,x,'linear',dxm(1)); % extrapolate to first dxm if out of bounds
if(My>1)
    %y = y./interp1(ym,dym,y,'linear',dym(1));
else
    %y = y/dym;
end
x = x/dxm(1);
y = y/dym(1);
xm = xm./dxm;
ym = ym./dym;
z = zeros(size(y));
e = zeros(size(y));
% compute basis function value at each output location and mult by amplitude
N = length(x);
for i=1:N
    % find basis functions that play role
    idx = find(abs(x(i)-xm)<=2);
    idy = find(abs(y(i)-ym)<=2);
    % saves computing alot of zeros
    if(length(idx)>0 & length(idy)>0)
        % evaluate at i
        fx = repmat(bspline_basis(x(i)-xm(idx))', length(idy),1);
        fy = repmat(bspline_basis(y(i)-ym(idy)) , 1,length(idx));
        %z(i) = col(am(idy,idx))'*f(:);
        %if(exist('aci'))
           %e(i) = ((col(aci(idy,idx))').^2)*(f(:).^2);
        %end
        % add bc contribution
        if(1)
            if(Mx>4 & My>4)
            %if(any(idx<3) | any(idy<3) | any(idx>(Mx-2)) | any(idy>(My-2)) & exist('bctype')) % first see if we need to bother
            if(any(idx<3) & exist('bctype') )
                % get x contribution
                fbx =  repmat(sparse(bspline_basis(x(i)-(xm(1)-1)))', length(idy),length(idx));
            elseif(any(idx>(Mx-2)) & exist('bctype'))
                fbx =  repmat(sparse(bspline_basis(x(i)-(xm(Mx)+1)))', length(idy),length(idx));
            else
                fbx = zeros(length(idy),length(idx));
            end
            if(any(idy<3) & exist('bctype'))
                % get y contribution
                fby =  repmat(sparse(bspline_basis(y(i)-(ym(1)-1))),length(idy),length(idx));
            elseif(any(idy>(My-2)) & exist('bctype'))
                fby =  repmat(sparse(bspline_basis(y(i)-(ym(My)+1))),length(idy),length(idx));
            else
                fby = zeros(length(idy),length(idx));
            end
            %abc = repmat(bcx(idx)',length(idy),1).*repmat(bcy(idy),1,length(idx)).*am(idy,idx);
            else
                fbx = zeros(length(idy),length(idx));
                fby = zeros(length(idy),length(idx));
            end
            
            try
                %z(i) = z(i) + abc(:)'*(fbx(:).*fby(:));
                f = (fx + repmat(bcx(idx)',length(idy),1).*fbx).*(fy + repmat(bcy(idy),1,length(idx)).*fby);
                z(i) = col(am(idy,idx))'*f(:);
                if(exist('aci'))
                    e(i) = (col(aci(idy,idx)).^2)'*(f(:).^2);
                end
            catch
                lasterr;
                keyboard
            end
            
            %if(exist('aci'))
                %abc = repmat(bcx(idx)',length(idy),1).*repmat(bcy(idy),1,length(idx)).*aci(idy,idx);
                %e(i) = e(i) + (abc(:).^2)'*[(fbx(:).*fby).^2];
            %end
        end
        
    else
        y(i) = 0;
    end
end


% need this
function c = col(c)
c = c(:);
return;