% test_bsplinecompute2d

dx = 1;
Nx = 21*1e1; % choose odd
dy = 2;
Ny = 41*1e1;
%Ny = 1;
x = [-0.5*(Nx-1):0.5*(Nx-1)]'*dx;
y = [-0.5*(Ny-1):0.5*(Ny-1)]'*dy;
[xo,yo] = meshgrid(x,y);
zo = exp(-(xo.^2)/(0.125*Nx^2) + -(yo.^2)/(0.125*Ny^2));%ones(Ny,Nx);%xo+yo; %rand(Nx*Ny,1);
%zo = 0*xo+1;
%zo = xo - (min(x));
%zo = yo - (min(y));
imagesc(zo); colorbar
w = 1;

dxm = 2*dx;
dym = 2*dy;
Nxm = 3+2*fix( (x(Nx)-x(1))/dxm/2 );
Nym = 3+2*fix( (y(Ny)-y(1))/dym/2 );
xm = [-0.5*(Nxm-1):0.5*(Nxm-1)]'*dxm;
ym = [-0.5*(Nym-1):0.5*(Nym-1)]'*dym;
if(0&Ny==1)
    Nym = 1;
    ym = 0;
end

lc = 0;
bctype = 0;
bctype =2; % 
bctype = -1; % do nothing

% works in 1-d,2d
%[am,aci,J] = bspline_compute2d_v03_a(xo(:),yo(:),zo(:),w,xm,ym,lc,bctype);

% working now
%[am,aci,J] = bspline_compute2d_v08(xo(:),yo(:),zo(:),w,xm,ym,lc,bctype);

% missing data
id = [1:5:(Ny*Nx)];
addpath('m:/plant/matlab/misc')
lc = 4;
bctype = 0;
[am,aci,J] = bspline_compute2d_v09(col(xo(id)),col(yo(id)),col(zo(id)),w,xm,ym,lc,bctype);

% show coeff
subplot 121
imagesc(xm,ym,am)
colorbar
subplot 122
imagesc(xm,ym,aci)
colorbar
hold on
%plot(xo(id),yo(id),'w.')
hold off


% check reconstruction
[zbs] = bspline_curve2d(xo(:),yo(:), xm,ym, am, bctype);
zbs = reshape(zbs, Ny,Nx);

subplot 121
imagesc(x,y,zo); caxis([-1 1]) 
hold on
%plot(xo(id),yo(id),'w.')
hold off
colorbar

subplot 122
imagesc(x,y,zbs)
caxis([-1 1])
colorbar
hold on
%plot(xo(id),yo(id),'w.')
hold off


% on 10 June 2005: reconstruction not quite right when bctype is used

%%%%%%%%%%%%%%%%%%%%%%%
% do it to some data
cd('E:\Users\Plant\Projects\02-interpolation\InterpNCEXBathy\ncexForD3D-7(newgrid28oct2004)')
load  outerswan_BSplineLocal.mat
w=0.1;
w = (1-NEi)/(1+w);
xm = Xlocal(1,1:2:end);
ym = Ylocal(1:2:end,1);
lc = 4;
bctype = 1;
[am,aci,J,b,r] = bspline_compute2d_v10(col(Xlocal),col(Ylocal),col(Zi),col(w),xm(:),ym(:),lc,bctype);

% show coeff
figure
subplot 121
imagesc(xm,ym,am)
colorbar
subplot 122
imagesc(xm,ym,aci)
colorbar
hold on
%plot(xo(id),yo(id),'w.')
hold off

% check reconstruction
[zbs] = bspline_curve2d(Xlocal(:),Ylocal(:), xm(:),ym(:), am, bctype);
zbs = reshape(zbs, size(Xlocal));

figure
x = Xlocal(1,:);
y = Ylocal(:,1);
subplot 121
imagesc(x,y,-Zbs); %caxis([-1 1]) 
hold on
%plot(xo(id),yo(id),'w.')
hold off
caxis([-500 10])
colorbar

subplot 122
imagesc(x,y,-zbs)
%caxis([-1 1])
caxis([-500 10])
colorbar
hold on
%plot(xo(id),yo(id),'w.')
hold off

figure
for i=1:length(y)
    plot(x, -Zi(i,:),'b.', x,-Zbs(i,:),'k', x, -zbs(i,:),'r', x, NEi(i,:)*10,'b--')
    grid on
    zoom on
    ylim([-500 10])
    pause
end
    














subplot 211
imagesc(am);
subplot 212
imagesc(aci)

% reconstruct
[zbs] = bspline_curve2d(x,y, xm,ym, reshape(am, My,Mx),bctype);
subplot 211
imagesc(reshape(z,21,11));
subplot 212
imagesc(reshape(zbs,21,11))


plot(1:21,reshape(z,21,11),'b',1:21,reshape(zbs,21,11),'r-' )


return

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
