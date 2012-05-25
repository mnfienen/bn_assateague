function Zi2 = splineFilter(Zi, Ei, dxm, dym)
% this function is a quick bspline filter to take
% a regular grid of depths (or other) and smooth it.
% dxm, dym are the subsample rates, e.g. dmx=2 estimates a spline parameter at every other grid point

splinedxm = [dxm, dym];
splinelc = 2; % smootheness control in regions with no data
splinebc = [1,1]; % 1st derivitive is zero on x boundaries

% guess at weights
EI_expected = (0.1).^2;
W = 1-(Ei)./(EI_expected+Ei);

ZB = mean(Zi(:));
ZP = Zi - ZB;

ZPS = bspline_pertgrid(ZP,W,splinebc,splinelc,splinedxm);

Zi2 = ZPS + ZB;
