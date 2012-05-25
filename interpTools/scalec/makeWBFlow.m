function wbflow = makeWBFlow(Ylocal, Nysmooth,mindyg)
wbflow = exp(-((Ylocal-Ylocal(1))/(1*Nysmooth*mindyg)).^2);
wbflow = wbflow + exp(-((Ylocal-Ylocal(end))/(1*Nysmooth*mindyg)).^2);
wbflow = 1-wbflow;
wbflow = wbflow/max(wbflow(:));

return


wb = 1-(Ei)./(wbfloor + Ei);
if(isflow>0)
    % change the spline type
    %splinebc = [1,10]; % different xshore,yshore bc
    % force to uniform
    %wb = wb.*(1-wbflow);
    wb = wb.*(wbflow);
end

[Zbs] = bspline_pertgrid(Zi.*(wbflow.^2),wb,splinebcpert,splinelc,splinedxm);
