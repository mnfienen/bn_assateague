% test_bsplinecompute2d

dx = 1;
Nx = 61*1e0; % choose odd
dy = 2;
Ny = 41*1e0;
%Ny = 1;
x = [-0.5*(Nx-1):0.5*(Nx-1)]'*dx;
y = [-0.5*(Ny-1):0.5*(Ny-1)]'*dy;
[xo,yo] = meshgrid(x,y);
%zo = exp(-(xo.^2)/(0.125*Nx^2) + -(yo.^2)/(0.125*Ny^2));%ones(Ny,Nx);%xo+yo; %rand(Nx*Ny,1);
%zo = 0*xo+1;
zo = xo - (mean(x));
%zo = yo - (min(y));
imagesc(zo); colorbar
w = 1;

dxm = 4*dx;
dym = 4*dy;
Nxm = 2+2*fix( (x(Nx)-x(1))/dxm/2 );
Nym = 2+2*fix( (y(Ny)-y(1))/dym/2 );
xm = [-0.5*(Nxm-1):0.5*(Nxm-1)]'*dxm;
ym = [-0.5*(Nym-1):0.5*(Nym-1)]'*dym;
if(0&Ny==1)
    Nym = 1;
    ym = 0;
end

% works in 1-d,2d
%[am,aci,J] = bspline_compute2d_v03_a(xo(:),yo(:),zo(:),w,xm,ym,lc,bctype);

% working now
%[am,aci,J] = bspline_compute2d_v08(xo(:),yo(:),zo(:),w,xm,ym,lc,bctype);

% missing data
addpath('m:/plant/matlab/misc')
id = [1:1:(Ny*Nx)];
lc = 0;
bctype =0;
[am,aci,J,b,r] = bspline_compute2d_v11b(col(xo(id)),col(yo(id)),col(zo(id)),w,xm,ym,lc,bctype);

% show coeff
subplot 121
imagesc(xm,ym,am)
colorbar
subplot 122
imagesc(xm,ym,aci)
colorbar
hold on
plot(xo(id),yo(id),'w.')
hold off

% show b
imagesc(xm,ym,reshape(b, Nym, Nxm))
hold on
plot(xo,yo,'w.')
hold off

% show slice
xslice = [(xm(1)-dxm):dx:(xm(end)+dxm)]';
yslice = 0*xslice;
[zbslice] = bspline_curve2d(xslice,yslice, xm,ym, am, bctype);
plot((xslice-xm(1))/dxm, zbslice, 'm', (xm-xm(1))/dxm, am(1+(Nym-1)/2,:),'mo', (x-xm(1))/dxm, zo(1,:),'b.'); grid on
grid on
hold on

% compare
lc = 0;
bctype =1;
xslice = [(xm(1)-dxm):dx:(xm(end)+dxm)]';
yslice = 0*xslice+x(1);
% 2d
[am,aci,J,b,r] = bspline_compute2d_v11b(col(xo(id)),col(yo(id)),col(zo(id)),w,xm,ym,lc,bctype);
[zbslice] = bspline_curve2d(xslice,yslice, xm,ym, am, bctype);

% 1d
[am1,aci1,J1,zm1,zs1] = bspline_compute(x(:),col(zo(1,:)),w,xm,dxm,lc,bctype);
[zbslice1] = bspline_curve(xslice, xm, am1, dxm, bctype);

%plot
%plot((x-xm(1))/dxm, zo(1,:),'b.', (xslice-xm(1))/dxm, zbslice1,'k', (xm-xm(1))/dxm, am1,'ko', (xslice-xm(1))/dxm,zbslice,'r', (xm-xm(1))/dxm,am(1+(Nym-1)/2,:),'ro'); grid on
plot((x-xm(1))/dxm, zo(1,:),'b.', (xslice-xm(1))/dxm, zbslice1,'k', (xm-xm(1))/dxm, am1,'ko', (xslice-xm(1))/dxm,zbslice,'r', (xm-xm(1))/dxm,am(1,:),'ro'); grid on




% check reconstruction
[zbs] = bspline_curve2d(xo(:),yo(:), xm,ym, am, bctype);
zbs = reshape(zbs, Ny,Nx);

subplot 121
imagesc(x,y,zo); 
caxis([min(zo(:)) max(zo(:))]) 
hold on
%plot(xo(id),yo(id),'w.')
hold off
colorbar

subplot 122
imagesc(x,y,zbs-zo)
caxis([min(zo(:)) max(zo(:))]) 
caxis([-1 1]) 
colorbar
hold on
%plot(xo(id),yo(id),'w.')
hold off


% on 10 June 2005: reconstruction not quite right when bctype is used

%%%%%%%%%%%%%%%%%%%%%%%
% do it to some data
addpath g:/matlab
pc_remotestartup
cd('E:\Users\Plant\Projects\02-interpolation\InterpNCEXBathy\ncexForD3D-7(newgrid28oct2004)')
load  outerswan_BSplineLocal.mat
w=0.5;
w = 1-NEi/(1+w);
dxm = 2*diff(Xlocal(1,1:2))
dym = 2*diff(Ylocal(1:2,1))
xm = (Xlocal(1)-dxm):dxm:(Xlocal(end)+dxm);
ym = (Ylocal(1)-dym):dym:(Ylocal(end)+dym);
lc = 4;
bctype = 1;
%
[am,aci,J,b,r] = bspline_compute2d_v11b(col(Xlocal),col(Ylocal),col(Zi),col(w),xm(:),ym(:),lc,bctype);
[zbs] = bspline_curve2d(Xlocal(:),Ylocal(:), xm(:),ym(:), am, bctype);
zbs = reshape(zbs, size(Xlocal));

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
    plot(x,-Zi(i,:),'b.',x,-Zi(i,:),'b-', x,-Zbs(i,:),'k', x, -zbs(i,:),'r', x, NEi(i,:)*10,'b--')
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
